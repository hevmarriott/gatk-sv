version 1.0

import "CollectCoverage.wdl" as cov
import "GATKSVGenotype.wdl" as svg
import "GermlineCNVTasks.wdl" as gcnv_tasks
import "GATKSVDepth.wdl" as gatksv_depth

workflow GATKSVJoinSamples {
  input {
    Array[File] vcfs
    String batch
    Array[String] samples
    Array[File] counts
    Array[File] cnmops_files

    File sr_file
    File pe_file
    File sample_mean_depth_file
    File sample_median_count_file
    File contig_list
    File ploidy_calls_tar

    # Filtering options
    File? inclusion_intervals_depth_only
    File? exclusion_intervals_depth_only
    Int min_size_depth_only = 5000
    Float min_overlap_fraction_depth_only = 0.5
    Boolean require_breakend_overlap_depth_only = false

    File? inclusion_intervals_non_depth_only
    File? exclusion_intervals_non_depth_only
    Int min_size_non_depth_only = 50
    Float min_overlap_fraction_non_depth_only = 0
    Boolean require_breakend_overlap_non_depth_only = true

    # Condense read counts
    Int small_cnv_condense_num_bins = 2
    Int small_cnv_condense_bin_size = 200
    Int large_cnv_condense_num_bins = 20
    Int large_cnv_condense_bin_size = 2000

    Int num_intervals_per_scatter = 10000

    Int large_cnv_padding = 1
    Int small_cnv_padding = 1
    Int depth_train_max_iter = 2000
    Int depth_predictive_samples = 100
    Int depth_predictive_iter = 10
    Int depth_discrete_samples = 1000

    Float? depth_mu_eps
    Float? depth_alpha_ref
    Float? depth_alpha_non_ref
    Float? depth_var_phi

    String depth_train_device = 'cpu'
    String depth_infer_device = 'cpu'

    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String linux_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_base_docker
    String gatk_docker
    String condense_counts_docker

    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_filter_depth
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_small_intervals
    RuntimeAttr? runtime_attr_intersect_intervals
    RuntimeAttr? runtime_attr_counts_to_intervals
    RuntimeAttr? runtime_attr_cluster
    RuntimeAttr? runtime_attr_posteriors
    RuntimeAttr? runtime_attr_condense_counts
    RuntimeAttr? runtime_attr_override_make_bincov
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_infer
    }

  scatter (i in range(length(vcfs))) {
    File vcf_indexes_ = vcfs[i] + ".tbi"
  }

  scatter (i in range(length(cnmops_files))) {
    File cnmops_file_indexes_ = cnmops_files[i] + ".tbi"
  }

  File sr_index_ = sr_file + ".tbi"
  File pe_index_ = pe_file + ".tbi"

  call MergeSVCalls {
    input:
      vcfs = vcfs,
      vcf_indexes = vcf_indexes_,
      cnmops_files = cnmops_files,
      cnmops_file_indexes = cnmops_file_indexes_,
      batch = batch,
      ref_fasta_dict = ref_fasta_dict,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_merge
  }

  call FilterVariants {
    input:
      vcf = MergeSVCalls.out,
      vcf_index = MergeSVCalls.out_index,
      inclusion_intervals_depth_only = inclusion_intervals_depth_only,
      exclusion_intervals_depth_only = exclusion_intervals_depth_only,
      min_size_depth_only = min_size_depth_only,
      min_overlap_fraction_depth_only = min_overlap_fraction_depth_only,
      require_breakend_overlap_depth_only = require_breakend_overlap_depth_only,
      min_size_non_depth_only = min_size_non_depth_only,
      min_overlap_fraction_non_depth_only = min_overlap_fraction_non_depth_only,
      require_breakend_overlap_non_depth_only = require_breakend_overlap_non_depth_only,
      inclusion_intervals_non_depth_only = inclusion_intervals_non_depth_only,
      exclusion_intervals_non_depth_only = exclusion_intervals_non_depth_only,
      basename = batch,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_filter_depth
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (contig in contigs) {
    call ClusterVariants {
      input:
        vcf = FilterVariants.out,
        vcf_index = FilterVariants.out_index,
        vid_prefix = "SV_" + contig + "_",
        contig = contig,
        sr_file = sr_file,
        sr_index = sr_index_,
        pe_file = pe_file,
        pe_index = pe_index_,
        sample_mean_depth_file = sample_mean_depth_file,
        ref_fasta_dict = ref_fasta_dict,
        batch = batch + "." + contig,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_cluster
    }
  }

  call svg.ConcatVcfs {
    input:
      vcfs = ClusterVariants.out,
      vcfs_idx = ClusterVariants.out_index,
      merge_sort = true,
      outfile_prefix = "~{batch}.clustered",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  call gatksv_depth.GATKSVDepth as DepthLarge {
    input:
      batch = batch,
      cnv_size_name = "large_cnv",
      vcf = ConcatVcfs.out,
      samples = samples,
      counts = counts,
      sample_median_count_file = sample_median_count_file,
      ploidy_calls_tar = ploidy_calls_tar,
      condense_num_bins = large_cnv_condense_num_bins,
      condense_bin_size = large_cnv_condense_bin_size,
      include_depth_only = true,
      cnv_size_conditional = ">= 5000",
      cnv_padding = large_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
      alpha_ref = depth_alpha_ref,
      train_max_iter = depth_train_max_iter,
      train_device = depth_train_device,
      predictive_samples = depth_predictive_samples,
      predictive_iter = depth_predictive_iter,
      discrete_samples = depth_discrete_samples,
      infer_device = depth_infer_device,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_fai = ref_fasta_fai,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      gatk_docker = gatk_docker,
      condense_counts_docker = condense_counts_docker,
      runtime_attr_small_intervals = runtime_attr_small_intervals,
      runtime_attr_override_make_bincov = runtime_attr_override_make_bincov,
      runtime_attr_intersect_intervals = runtime_attr_intersect_intervals,
      runtime_attr_counts_to_intervals = runtime_attr_counts_to_intervals,
      runtime_attr_condense_counts = runtime_attr_condense_counts,
      runtime_attr_scatter = runtime_attr_scatter,
      runtime_attr_train = runtime_attr_train,
      runtime_attr_infer = runtime_attr_infer,
      runtime_attr_concat = runtime_attr_concat
  }

  call gatksv_depth.GATKSVDepth as DepthSmall {
    input:
      batch = batch,
      cnv_size_name = "small_cnv",
      vcf = ConcatVcfs.out,
      samples = samples,
      counts = counts,
      sample_median_count_file = sample_median_count_file,
      ploidy_calls_tar = ploidy_calls_tar,
      condense_num_bins = small_cnv_condense_num_bins,
      condense_bin_size = small_cnv_condense_bin_size,
      include_depth_only = false,
      cnv_size_conditional = "< 5000",
      cnv_padding = small_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
      mu_eps = depth_mu_eps,
      alpha_ref = depth_alpha_ref,
      alpha_non_ref = depth_alpha_non_ref,
      var_phi = depth_var_phi,
      train_max_iter = depth_train_max_iter,
      train_device = depth_train_device,
      predictive_samples = depth_predictive_samples,
      predictive_iter = depth_predictive_iter,
      discrete_samples = depth_discrete_samples,
      infer_device = depth_infer_device,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_fai = ref_fasta_fai,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      gatk_docker = gatk_docker,
      condense_counts_docker = condense_counts_docker,
      runtime_attr_small_intervals = runtime_attr_small_intervals,
      runtime_attr_override_make_bincov = runtime_attr_override_make_bincov,
      runtime_attr_intersect_intervals = runtime_attr_intersect_intervals,
      runtime_attr_counts_to_intervals = runtime_attr_counts_to_intervals,
      runtime_attr_condense_counts = runtime_attr_condense_counts,
      runtime_attr_scatter = runtime_attr_scatter,
      runtime_attr_train = runtime_attr_train,
      runtime_attr_infer = runtime_attr_infer,
      runtime_attr_concat = runtime_attr_concat
  }

  call AggregateDepth {
    input:
      vcf = ConcatVcfs.out,
      vcf_index = ConcatVcfs.out_index,
      depth_posterior_vcfs = [DepthLarge.out, DepthSmall.out],
      depth_posterior_vcfs_indexes = [DepthLarge.out_index, DepthSmall.out_index],
      ploidy_calls_tar = ploidy_calls_tar,
      ref_fasta_dict = ref_fasta_dict,
      batch = batch,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_posteriors
  }

  output {
    File joined_vcf = AggregateDepth.out
    File joined_vcf_index = AggregateDepth.out_index

    File large_cnv_depth_vcf = DepthLarge.out
    File large_cnv_depth_vcf_index = DepthLarge.out_index
    File large_cnv_depth_counts = DepthLarge.depth_file
    File large_cnv_depth_counts_index = DepthLarge.depth_file_index

    File small_cnv_depth_vcf = DepthSmall.out
    File small_cnv_depth_vcf_index = DepthSmall.out_index
    File small_cnv_depth_counts = DepthSmall.depth_file
    File small_cnv_depth_counts_index = DepthSmall.depth_file_index
  }
}

task AggregateDepth {
  input {
    File vcf
    File vcf_index
    Array[File] depth_posterior_vcfs
    Array[File] depth_posterior_vcfs_indexes
    File ploidy_calls_tar
    File ref_fasta_dict
    String batch
    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{batch}.aggregated.vcf.gz"
    File out_index = "~{batch}.aggregated.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    mkdir ploidy-calls
    tar xzf ~{ploidy_calls_tar} -C ploidy-calls
    ls ploidy-calls/SAMPLE_*/contig_ploidy.tsv > ploidy_files.list

    # Create arguments file
    echo "--cnv-intervals-vcf ~{sep=" --cnv-intervals-vcf " depth_posterior_vcfs}" > args.txt
    while read line; do
      echo "--ploidy-calls-file $line" >> args.txt
    done < ploidy_files.list

    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVAggregateDepth \
      --arguments_file args.txt \
      --variant ~{vcf} \
      --output ~{batch}.aggregated.vcf.gz \
      --sequence-dictionary ~{ref_fasta_dict} \
      --genotype-depth-calls

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ClusterVariants {
  input {
    File vcf
    File vcf_index
    String vid_prefix
    String contig
    File sr_file
    File sr_index
    File pe_file
    File pe_index
    File sample_mean_depth_file
    File ref_fasta_dict
    String batch
    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: {
      localization_optional: true
    }
    sr_file: {
      localization_optional: true
    }
    pe_file: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 7.5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{batch}.clustered.vcf.gz"
    File out_index = "~{batch}.clustered.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVCluster \
      -V ~{vcf} \
      -O ~{batch}.clustered.vcf.gz \
      -L ~{contig} \
      --sequence-dictionary ~{ref_fasta_dict} \
      --split-reads-file ~{sr_file} \
      --discordant-pairs-file ~{pe_file} \
      --sample-coverage ~{sample_mean_depth_file} \
      --variant-prefix ~{vid_prefix}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task MergeSVCalls {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[File] cnmops_files
    Array[File] cnmops_file_indexes
    File ref_fasta_dict
    String batch
    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 7.5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{batch}.merged.vcf.gz"
    File out_index = "~{batch}.merged.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Create arguments file
    touch args.txt
    while read line; do
      echo "--cnmops $line" >> args.txt
    done < ~{write_lines(cnmops_files)}

    while read line; do
      echo "-V $line" >> args.txt
    done < ~{write_lines(vcfs)}

    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" MergeSVCalls \
      --arguments_file args.txt \
      --sequence-dictionary ~{ref_fasta_dict} \
      --output ~{batch}.merged.vcf.gz \
      --ignore-dict
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task MergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_name}"
    File out_index = "~{output_name}.tbi"
  }
  command <<<

    set -euo pipefail
    bcftools merge --file-list ~{write_lines(vcfs)} -o ~{output_name} --output-type z
    tabix ~{output_name}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ShardVcf {
  input {
    File pesr_vcf
    Int records_per_shard
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Array[String] svtypes = ["DEL", "DUP", "INV", "INS"]
  Array[String] bnd_strands = ["++", "+-", "-+", "--"]
  Array[String] bnd_strand_labels = ["plus_plus", "plus_minus", "minus_plus", "minus_minus"]

  output {
    Array[File] out = glob("*.shard_*.vcf.gz")
  }
  command <<<

    set -euo pipefail
    sgrep() { grep "$@" || test $? = 1; }
    NAME=$(basename ~{pesr_vcf} .vcf.gz)
    SVTYPES=(~{sep=" " svtypes})
    zcat ~{pesr_vcf} | grep ^# > header
    for svtype in ${SVTYPES[@]}; do
      zcat ~{pesr_vcf} | grep -v ^# | sgrep -w "<${svtype}>" > records
      NUM_RECORDS=$(cat records | wc -l)
      if [ "${NUM_RECORDS}" -gt "0" ]; then
        CHUNKS=$(python -c "from math import ceil; print(ceil($NUM_RECORDS/float(~{records_per_shard})))")
        split -a9 -d -n l/$CHUNKS records records.shard_
        i=0
        for chunk in records.shard_*; do
          padded=`printf %09d $i`
          cat header $chunk | bgzip -c > $NAME.${svtype}.shard_${padded}.vcf.gz
          i=$((i+1))
        done
        rm records records.shard_*
      else
        echo "No records of type ${svtype} found"
      fi
    done

    # Handle BNDs separately
    STRANDS=(~{sep=" " bnd_strands})
    STRAND_STR=(~{sep=" " bnd_strand_labels})
    svtype="BND"
    num_strand_types=${#STRANDS[*]}
    for (( j=0; j<=$(( $num_strand_types -1 )); j++ )); do
      strand=${STRANDS[$j]}
      strand_str=${STRAND_STR[$j]}
      zcat ~{pesr_vcf} | grep -v ^# | sgrep -w "<${svtype}>" | sgrep "STRANDS=${strand}" > records
      NUM_RECORDS=$(cat records | wc -l)
      if [ "${NUM_RECORDS}" -gt "0" ]; then
        CHUNKS=$(python -c "from math import ceil; print(ceil($NUM_RECORDS/float(~{records_per_shard})))")
        split -a9 -d -n l/$CHUNKS records records.shard_
        i=0
        for chunk in records.shard_*; do
          padded=`printf %09d $i`
          cat header $chunk | bgzip -c > $NAME.${svtype}.${strand_str}.shard_${padded}.vcf.gz
          i=$((i+1))
        done
        rm records records.shard_*
      else
        echo "No records of type ${svtype} found"
      fi
    done

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task FilterVariants {
  input {
    File vcf
    File vcf_index

    File? inclusion_intervals_depth_only
    File? exclusion_intervals_depth_only
    Int min_size_depth_only
    Float min_overlap_fraction_depth_only
    Boolean require_breakend_overlap_depth_only

    File? inclusion_intervals_non_depth_only
    File? exclusion_intervals_non_depth_only
    Int min_size_non_depth_only
    Float min_overlap_fraction_non_depth_only
    Boolean require_breakend_overlap_non_depth_only

    String basename
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{basename}.filtered.vcf.gz"
    File out_index = "~{basename}.filtered.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    gatk --java-options -Xmx~{java_mem_mb}M SVSelectVariants \
      -V ~{vcf} \
      -O ~{basename}.depth_filtered.vcf.gz \
      --depth-only \
      ~{"-L " + inclusion_intervals_depth_only} \
      ~{"-XL " + exclusion_intervals_depth_only} \
      --min-size ~{min_size_depth_only} \
      --min-overlap-fraction ~{min_overlap_fraction_depth_only} \
      ~{if require_breakend_overlap_depth_only then "--require-breakend-overlap" else ""}

    gatk --java-options -Xmx~{java_mem_mb}M SVSelectVariants \
      -V ~{vcf} \
      -O ~{basename}.non_depth_filtered.vcf.gz \
      --non-depth-only \
      ~{"-L " + inclusion_intervals_non_depth_only} \
      ~{"-XL " + exclusion_intervals_non_depth_only} \
      --min-size ~{min_size_non_depth_only} \
      --min-overlap-fraction ~{min_overlap_fraction_non_depth_only} \
      ~{if require_breakend_overlap_non_depth_only then "--require-breakend-overlap" else ""}

    bcftools concat \
      --allow-overlaps \
      ~{basename}.depth_filtered.vcf.gz ~{basename}.non_depth_filtered.vcf.gz \
      | bgzip > ~{basename}.filtered.vcf.gz
    tabix ~{basename}.filtered.vcf.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
