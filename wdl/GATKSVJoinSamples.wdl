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
    Array[File] large_gcnv_interval_vcfs

    File sr_file
    File pe_file
    File sample_coverage_file
    File gatk_sv_cluster_exclude_intervals
    File exclude_intervals
    File contig_list
    File ploidy_calls_tar

    Int min_depth_only_size = 5000

    # Condense read counts
    Int small_cnv_condense_num_bins = 2
    Int small_cnv_condense_bin_size = 200
    Int large_cnv_condense_num_bins = 20
    Int large_cnv_condense_bin_size = 2000

    Int num_intervals_per_scatter = 10000

    Int large_cnv_padding = 100000
    Int small_cnv_padding = 1000
    Int depth_train_max_iter = 2000
    Int depth_predictive_samples = 100
    Int depth_predictive_iter = 10
    Int depth_discrete_samples = 1000

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

  scatter (i in range(length(large_gcnv_interval_vcfs))) {
    File large_gcnv_interval_vcf_indexes_ = large_gcnv_interval_vcfs[i] + ".tbi"
  }

  File sr_index_ = sr_file + ".tbi"
  File pe_index_ = pe_file + ".tbi"
  File exclude_intervals_index_ = exclude_intervals + ".tbi"

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

  Array[String] contigs = read_lines(contig_list)

  scatter (contig in contigs) {
    call ClusterVariants {
      input:
        vcf = MergeSVCalls.out,
        vcf_index = MergeSVCalls.out_index,
        vid_prefix = "SV_" + contig + "_",
        contig = contig,
        sr_file = sr_file,
        sr_index = sr_index_,
        pe_file = pe_file,
        pe_index = pe_index_,
        sample_coverage_file = sample_coverage_file,
        gatk_sv_cluster_exclude_intervals = gatk_sv_cluster_exclude_intervals,
        exclude_intervals = exclude_intervals,
        exclude_intervals_index = exclude_intervals_index_,
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

  call FilterDepthOnlyBySize {
    input:
      vcf = ConcatVcfs.out,
      vcf_index = ConcatVcfs.out_index,
      size = min_depth_only_size,
      basename = batch,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_filter_depth
  }

  call gatksv_depth.GATKSVDepth as DepthLarge {
    input:
      batch = batch,
      cnv_size_name = "large_cnv",
      vcf = FilterDepthOnlyBySize.out,
      samples = samples,
      counts = counts,
      sample_coverage_file = sample_coverage_file,
      ploidy_calls_tar = ploidy_calls_tar,
      condense_num_bins = large_cnv_condense_num_bins,
      condense_bin_size = large_cnv_condense_bin_size,
      include_depth_only = true,
      cnv_size_conditional = ">= 5000",
      cnv_padding = large_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
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
      vcf = FilterDepthOnlyBySize.out,
      samples = samples,
      counts = counts,
      sample_coverage_file = sample_coverage_file,
      ploidy_calls_tar = ploidy_calls_tar,
      condense_num_bins = small_cnv_condense_num_bins,
      condense_bin_size = small_cnv_condense_bin_size,
      include_depth_only = false,
      cnv_size_conditional = "< 5000",
      cnv_padding = small_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
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

  call CopyNumberPosteriors {
    input:
      vcf = FilterDepthOnlyBySize.out,
      vcf_index = FilterDepthOnlyBySize.out_index,
      gcnv_intervals_vcfs = [DepthLarge.out, DepthLarge.out],
      gcnv_intervals_vcf_indexes = [DepthLarge.out_index, DepthLarge.out_index],
      ploidy_calls_tar = ploidy_calls_tar,
      ref_fasta_dict = ref_fasta_dict,
      batch = batch,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_posteriors
  }

  output {
    File out = CopyNumberPosteriors.out
    File out_index = CopyNumberPosteriors.out_index
  }
}

task CopyNumberPosteriors {
  input {
    File vcf
    File vcf_index
    Array[File] gcnv_intervals_vcfs
    Array[File] gcnv_intervals_vcf_indexes
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
    echo "--cnv-intervals-vcf ~{sep=" --cnv-intervals-vcf " gcnv_intervals_vcfs}" > args.txt
    while read line; do
      echo "--ploidy-calls-file $line" >> args.txt
    done < ploidy_files.list

    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVCopyNumberPosteriors \
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
    File sample_coverage_file
    File gatk_sv_cluster_exclude_intervals
    File exclude_intervals
    File exclude_intervals_index
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
    disk_gb: 100,
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
      -XL ~{gatk_sv_cluster_exclude_intervals} \
      -XL ~{exclude_intervals} \
      --sequence-dictionary ~{ref_fasta_dict} \
      --split-reads-file ~{sr_file} \
      --discordant-pairs-file ~{pe_file} \
      --sample-coverage ~{sample_coverage_file} \
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

task SVTrainDepth {
  input {
    File depth_file
    File intervals
    File ploidy_calls_tar
    File sample_coverage_file
    File ref_dict
    String model_name
    String gatk_docker
    String device
    Int? max_iter
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    depth_file: {
           localization_optional: true
         }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 15,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{model_name}.cnv_model.tar.gz"
  }
  command <<<
    set -euo pipefail

    mkdir ploidy-calls
    tar xzf ~{ploidy_calls_tar} -C ploidy-calls
    ls ploidy-calls/SAMPLE_*/contig_ploidy.tsv > ploidy_files.list

    # Create arguments file
    while read line; do
    echo "--ploidy-calls-file $line" >> args.txt
    done < ploidy_files.list

    mkdir svmodel
    gatk --java-options -Xmx~{java_mem_mb}M SVTrainDepth \
      -L ~{intervals} \
      --depth-file ~{depth_file} \
      --coverage-file ~{sample_coverage_file} \
      --output-name ~{model_name} \
      --output-dir svmodel \
      --arguments_file args.txt \
      --sequence-dictionary ~{ref_dict} \
      --jit \
      ~{"--max-iter " + max_iter}

    tar czf ~{model_name}.cnv_model.tar.gz svmodel/*
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

task SVInferDepth {
  input {
    File model_tar
    File ref_dict
    Int predictive_samples
    Int predictive_iter
    Int discrete_samples
    String model_name
    String output_vcf_filename
    String gatk_docker
    String device
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: 10,
                               boot_disk_gb: 15,
                               preemptible_tries: 3,
                               max_retries: 0
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_vcf_filename}"
    File out_index = "~{output_vcf_filename}.tbi"
  }
  command <<<

    set -eo pipefail
    mkdir svmodel
    tar xzf ~{model_tar} svmodel/

    gatk --java-options -Xmx~{java_mem_mb}M SVInferDepth \
      --output ~{output_vcf_filename} \
      --predictive-samples ~{predictive_samples} \
      --predictive-iter ~{predictive_iter} \
      --discrete-samples ~{discrete_samples} \
      --model-name ~{model_name} \
      --model-dir svmodel \
      --sequence-dictionary ~{ref_dict} \
      --jit

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

task SmallCNVIntervals {
  input {
    File vcf
    File vcf_index
    File ref_fasta_fai
    Int size = 5000
    Int padding = 1000
    String batch
    String sv_pipeline_docker
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

  output {
    File out = "~{batch}.small_cnv_intervals.bed"
  }
  command <<<
    set -euo pipefail
    python > small_cnvs.bed <<EOF
import sys
from pysam import VariantFile
vcf = VariantFile('~{vcf}')
types = set(['DEL', 'DUP', 'BND'])
for record in vcf.fetch():
  if record.info['SVTYPE'] in types \
    and record.info['ALGORITHMS'] != 'depth' \
    and record.chrom == record.info['CHR2'] \
    and record.info['STRANDS'][0] != record.info['STRANDS'][1] \
    and record.stop - record.pos < ~{size}:
    fields = [record.chrom, str(record.pos), str(record.stop)]
    print('\t'.join(fields))
EOF

    bedtools slop -b ~{padding} -i small_cnvs.bed -g ~{ref_fasta_fai} \
      | bedtools merge \
      > ~{batch}.small_cnv_intervals.bed

  >>>
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

task IntersectCountsWithIntervals {
  input {
    File counts
    File interval_list
    String output_name
    Boolean gzip = true
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
    File out = "~{output_name}.gz"
  }
  command <<<

    set -euo pipefail
    zgrep -B9999999999 -m1 -v "^@" ~{counts} > ~{output_name}
    zgrep -v "^@" ~{counts} | tail -n +2 | bedtools intersect -wa -sorted -u -a stdin -b ~{interval_list} >> ~{output_name}
    if ~{gzip}; then
      bgzip ~{output_name}
    fi

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


task FilterDepthOnlyBySize {
  input {
    File vcf
    File vcf_index
    Int size
    String basename
    String gatk_docker
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{basename}.depth_filtered.vcf.gz"
    File out_index = "~{basename}.depth_filtered.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    gatk --java-options -Xmx~{java_mem_mb}M SelectVariants \
      -V ~{vcf} \
      -O ~{basename}.depth_filtered.vcf.gz \
      -select "ALGORITHMS == 'depth' && SVLEN <= ~{size}" \
      --invertSelect

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