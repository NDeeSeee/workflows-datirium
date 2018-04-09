cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement


inputs:
  upstream_fastq:
    type: File
  downstream_fastq:
    type: File?
  rsem_indices_folder:
    type: Directory
  chr_length_file:
    type: File
  paired:
    type: boolean?
  threads:
    type: int?


outputs:

  rsem_isoforms:
    type: File
    outputSource: rsem_calculate_expression/isoform_results_file
  rsem_genes:
    type: File
    outputSource: rsem_calculate_expression/gene_results_file
  bam_file:
    type: File
    outputSource: rsem_calculate_expression/genome_sorted_bam_bai_pair
  bamtools_log:
    type: File
    outputSource: bamtools_stats/log_file
  bed:
    type: File
    outputSource: bam_to_bigwig/bedgraph_file
  bigwig:
    type: File
    outputSource: bam_to_bigwig/bigwig_file


steps:

  rsem_calculate_expression:
    run: ../tools/rsem-calculate-expression.cwl
    in:
      upstream_read_file: upstream_fastq
      downstream_read_file: downstream_fastq
      indices_folder: rsem_indices_folder
      bowtie2:
        default: true
      sort_bam_by_coordinate:
        default: true
      output_genome_bam:
        default: true
      threads: threads
    out:
      - isoform_results_file
      - gene_results_file
      - genome_sorted_bam_bai_pair

  bamtools_stats:
    run: ../tools/bamtools-stats.cwl
    in:
      bam_file: rsem_calculate_expression/genome_sorted_bam_bai_pair
    out:
     - mapped_reads_number
     - log_file

  bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: rsem_calculate_expression/genome_sorted_bam_bai_pair
      chrom_length_file: chr_length_file
      mapped_reads_number: bamtools_stats/mapped_reads_number
      pairchip: paired
    out:
      - bigwig_file
      - bedgraph_file