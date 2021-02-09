version 1.0

import "Structs.wdl"
import "Tasks0506.wdl" as t0506

workflow GATKSVPreprocessSample {
  input {
    # Sample data (must provide at least one file)
    String sample_id
    File? manta_vcf
    File? melt_vcf
    File? wham_vcf
    Array[File] cnv_beds
    # Note this will not work on gCNV VCFs generated with GATK < v4.1.5.0 ! Use cnv_beds instead
    File? gcnv_segments_vcf

    # Filtering options
    Int min_svsize = 50
    Int gcnv_min_qs = 80

    # Reference index
    File ref_fasta_fai

    # Docker
    String sv_pipeline_docker

    # VM resource options
    RuntimeAttr? runtime_attr_override
  }

  Array[File?] optional_sv_vcfs_ = [manta_vcf, melt_vcf, wham_vcf]
  Array[String] optional_sv_algorithms_ = ["manta", "melt", "wham"]
  scatter (i in range(length(optional_sv_vcfs_))) {
    if (defined(optional_sv_vcfs_[i])) {
      File scattered_sv_vcfs_ = select_first([optional_sv_vcfs_[i]])
      String scattered_sv_algorithms_ = optional_sv_algorithms_[i]
    }
  }
  Array[File] sv_vcfs_ = select_all(scattered_sv_vcfs_)
  Array[String] sv_algorithms_ = select_all(scattered_sv_algorithms_)

  call StandardizeVcfs {
    input:
      vcfs = sv_vcfs_,
      algorithms = sv_algorithms_,
      sample = sample_id,
      gcnv_segments_vcf = gcnv_segments_vcf,
      cnv_beds = cnv_beds,
      output_basename = "~{sample_id}.preprocess_sample",
      gcnv_min_qs = gcnv_min_qs,
      min_svsize = min_svsize,
      ref_fasta_fai = ref_fasta_fai,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File out = StandardizeVcfs.out
    File out_index = StandardizeVcfs.out_index
  }
}

task StandardizeVcfs {
  input {
    Array[File] vcfs
    Array[String] algorithms
    File? gcnv_segments_vcf
    Array[File] cnv_beds
    String sample
    String output_basename
    Int gcnv_min_qs
    Int min_svsize
    File ref_fasta_fai
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_vcfs = length(vcfs)

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 3
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euxo pipefail

    ############################################################
    # Filter and standardize gCNV segments VCF
    ############################################################

    if ~{defined(gcnv_segments_vcf)}; then
      # Force GT type to String to avoid a bcftools bug
      tabix ~{gcnv_segments_vcf}
      zcat ~{gcnv_segments_vcf} \
        | sed 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' \
        | bgzip \
        > gcnv_reheadered.vcf.gz
      tabix gcnv_reheadered.vcf.gz

      # With older gCNV versions, piping vcf directly into bcftools gives this error:
      # [W::vcf_parse] Contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)
      # Note use of --no-version to avoid header timestamp, which breaks call caching
      bcftools view \
        --no-version \
        -s ~{sample} \
        -i 'ALT!="." && QS>~{gcnv_min_qs}' \
        -O z \
        -o gcnv_filtered.vcf.gz \
        gcnv_reheadered.vcf.gz
      tabix gcnv_filtered.vcf.gz

      # Standardize by adding required INFO fields
      # Note this will not work on VCFs generated with GATK < v4.1.5.0 !
      bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%ALT\n' gcnv_filtered.vcf.gz \
        | awk -F "\t" -v OFS="\t" '{
          if ($5=="<DEL>")  {
            svtype="DEL"; strands="+-";
          } else if ($5=="<DUP>") {
            svtype="DUP"; strands="-+";
          } else {
            svtype="."; strands=".";
          }
          print $1,$2,$3,$4,$3-$2+1,"depth",svtype,strands
        }' \
        | bgzip \
        > ann.bed.gz
      tabix -p bed ann.bed.gz

      echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">' > header_lines.txt
      echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> header_lines.txt
      echo '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness [++,+-,-+,--]">' >> header_lines.txt
      echo '##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">' >> header_lines.txt

      bcftools annotate \
        --no-version \
        -a ann.bed.gz \
        -c CHROM,POS,END,ID,SVLEN,ALGORITHMS,SVTYPE,STRANDS \
        -h header_lines.txt \
        -O z \
        -o gcnv.vcf.gz \
        gcnv_filtered.vcf.gz
      tabix gcnv.vcf.gz
      echo "gcnv.vcf.gz" > vcfs.list
    fi

    ############################################################
    # Standardize SV caller VCFs
    ############################################################

    vcfs=(~{sep=" " vcfs})
    algorithms=(~{sep=" " algorithms})
    for (( i=0; i<~{num_vcfs}; i++ )); do
      vcf=${vcfs[$i]}
      algorithm=${algorithms[$i]}
      tabix $vcf
      svtk standardize \
        --sample-names ~{sample} \
        --prefix ${algorithm}_~{sample} \
        --contigs ~{ref_fasta_fai} \
        --min-size ~{min_svsize} \
        $vcf unsorted.vcf ${algorithm}
      vcf-sort -c unsorted.vcf | bgzip > sorted.${algorithm}.vcf.gz
      tabix sorted.${algorithm}.vcf.gz
      echo "sorted.${algorithm}.vcf.gz" >> vcfs.list
    done

    ############################################################
    # Convert bed files into a standardized VCF
    ############################################################

    beds=(~{sep=" " cnv_beds})
    if ~{length(cnv_beds) > 0}; then
      echo "~{sample}" > samples.list
      # filter, concat, and header a combined bed
      zcat ~{sep=" " cnv_beds} \
        | awk -F "\t" -v OFS="\t" '{if ($5=="~{sample}") print}' \
        | grep -v ^# \
        | sort -k1,1V -k2,2n \
        > cnvs.bed
      # head command triggers pipefail
      set +o pipefail
      zcat ${beds[0]} | head -n1 > header.txt
      set -o pipefail
      cat header.txt cnvs.bed > cnvs.headered.bed

      # note svtk generates an index automatically
      svtk rdtest2vcf --contigs ~{ref_fasta_fai} cnvs.headered.bed samples.list cnvs.vcf.gz
      echo "cnvs.vcf.gz" >> vcfs.list
    fi

    ############################################################
    # Combine VCFs
    ############################################################

    bcftools concat -a --output-type z --file-list vcfs.list --output ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  >>>

  output {
    File out = "~{output_basename}.vcf.gz"
    File out_index = "~{output_basename}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
