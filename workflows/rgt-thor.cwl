cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  first_biological_condition:
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-pe.cwl"
  second_biological_condition:
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-pe.cwl"
  genome_indices:
  - "genome-indices.cwl"
  - "https://github.com/datirium/workflows/workflows/genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  bambai_pair_cond_1:
    type: File[]
    secondaryFiles:
      - .bai
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 1"
    doc: "Coordinate sorted BAM alignment and index BAI files for the first biological condition"
    'sd:upstreamSource': "first_biological_condition/bambai_pair"
    'sd:localLabel': true

  bambai_pair_cond_2:
    type: File[]
    secondaryFiles:
      - .bai
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 2"
    doc: "Coordinate sorted BAM alignment and index BAI files for the second biological condition"
    'sd:upstreamSource': "second_biological_condition/bambai_pair"
    'sd:localLabel': true

  alias_cond_1:
    type: string?
    default: "condition_1"
    label: "Name for condition 1"
    doc: "Name to be displayed for condition 1"
    'sd:layout':
      advanced: true

  alias_cond_2:
    type: string?
    default: "condition_2"
    label: "Name for condition 2"
    doc: "Name to be displayed for condition 2"
    'sd:layout':
      advanced: true

  chrom_length_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Chromosome length file"
    doc: "Chromosome length file"
    'sd:upstreamSource': "genome_indices/chrom_length"

  annotation_file:
    type: File
    label: "Genome annotation"
    format: "http://edamontology.org/format_3475"
    doc: "Genome annotation file in TSV format"
    'sd:upstreamSource': "genome_indices/annotation"

  merge_peaks:
    type: boolean?
    default: true
    label: "Merge peaks closer than fragment size"
    doc: "Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data)"
    'sd:layout':
      advanced: true

  remove_duplicates:
    type: boolean?
    default: false
    label: "Remove the duplicate reads"
    doc: "Remove the duplicate reads"
    'sd:layout':
      advanced: true

  housekeeping_genes_bed_file:
    type: File?
    format: "http://edamontology.org/format_3003"
    label: "Housekeeping genes file"
    doc: "Define housekeeping genes (BED format) used for normalizing"
    'sd:layout':
     advanced: true

  deadzones_bed_file:
    type: File?
    format: "http://edamontology.org/format_3003"
    label: "Dead zones file"
    doc: "Define blacklisted genomic regions avoided for analysis"
    'sd:layout':
     advanced: true

  pvalue_cutoff:
    type: float?
    default: 0.1
    label: "P-value cutoff for peak detection"
    doc: "P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: 0.1]"
    'sd:layout':
      advanced: true

  bin_size:
    type: int?
    default: 100
    label: "Size of underlying bins for creating the signal"
    doc: "Size of underlying bins for creating the signal"
    'sd:layout':
      advanced: true

  # no_correction:
  #   type: boolean?
  #   default: false
  #   label: "Skip p-value correction"
  #   doc: "Do not use multipe test correction for p-values (Benjamini/Hochberg)"
  #   'sd:layout':
  #     advanced: true

  extension_size:
    type:
      - "null"
      - string
      - int[]
    label: "Comma-separated list of read extension sizes (provide value for every sample)"
    doc: |
      Read's extension size for BAM files (comma separated list for each BAM file in config file).
      If option is not chosen, estimate extension sizes
    'sd:layout':
      advanced: true


outputs:

  diffpeaks_bed_file:
    type: File
    format: "http://edamontology.org/format_3004"
    label: "Estimated differential peaks"
    doc: "Estimated differential peaks, bigBed"
    outputSource: bed_to_bigbed/bigbed_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        format: 'bigbed'
        name: "Differential peaks"
        height: 120

  diffpeaks_annotated_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Estimated differential peaks with assigned genes"
    doc: "File contains nearest gene information for the differential peaks BED file generated by rgt-THOR"
    outputSource: restore_columns/output_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differential Peak Calling'
        Title: 'Differential Peaks rgt-THOR Results'

  cond_1_bigwig_file:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "First biological condition ChIP-seq signals"
    doc: "Postprocessed ChIP-seq signals from the first biological condition samples"
    outputSource: thor/cond_1_bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Biological condition 1"
        height: 120

  cond_2_bigwig_file:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Second biological condition ChIP-seq signals"
    doc: "Postprocessed ChIP-seq signals from the second biological condition samples"
    outputSource: thor/cond_2_bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Biological condition 2"
        height: 120

  thor_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "rgt-THOR stderr log"
    doc: "rgt-THOR stderr log"
    outputSource: thor/stderr_log


steps:

  thor:
    run: ../tools/rgt-thor.cwl
    in:
      bambai_pair_cond_1: bambai_pair_cond_1
      bambai_pair_cond_2: bambai_pair_cond_2
      chrom_length_file: chrom_length_file
      merge_peaks: merge_peaks
      housekeeping_genes_bed_file: housekeeping_genes_bed_file
      deadzones_bed_file: deadzones_bed_file
      pvalue_cutoff: pvalue_cutoff
      extension_size: extension_size
      # no_correction: no_correction
      remove_duplicates: remove_duplicates
      bin_size: bin_size
    out:
      - diffpeaks_bed_file
      - cond_1_bigwig_file
      - cond_2_bigwig_file
      - stderr_log

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: thor/diffpeaks_bed_file
      script:
        default: >
          cat $0 | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out: [output_file]

  assign_genes:
      run: ../tools/iaintersect.cwl
      in:
        input_filename: filter_columns/output_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file]

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: [assign_genes/result_file, thor/diffpeaks_bed_file]
      param: [alias_cond_1, alias_cond_2]
      script:
        default: |
          NAME_1=$2
          NAME_2=$3
          cat $0 | grep -v "start" | sort -k 11n > sorted_iaintersect_result.tsv
          cat $1 | tr ";" "\t" > thor_result.tsv
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tchrom\tstart\tend\tlength\tregion\tname\tscore\tcondition\tcolor\t${NAME_1}_counts\t${NAME_2}_counts\t-log10(pvalue)" > `basename $0`;
          cat sorted_iaintersect_result.tsv | paste - thor_result.tsv | cut -f 1-9,15,19-21,24,26-28 >> `basename $0`
          rm sorted_iaintersect_result.tsv thor_result.tsv
    out: [output_file]

  sort_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: thor/diffpeaks_bed_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: sort_bed/sorted_file
      bed_type:
        default: "bed4+7"
      chrom_length_file: chrom_length_file
      output_filename:
        source: sort_bed/sorted_file
        valueFrom: $(self.basename + ".bigBed")
    out: [bigbed_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "THOR - differential peak calling of ChIP-seq signals with replicates"
label: "THOR - differential peak calling of ChIP-seq signals with replicates"
s:alternateName: "THOR is an HMM-based approach to detect and analyze differential peaks in two sets of ChIP-seq data from distinct biological conditions with replicates"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/rgt-thor.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
  - class: s:Organization
    s:legalName: "Cincinnati Children's Hospital Medical Center"
    s:location:
    - class: s:PostalAddress
      s:addressCountry: "USA"
      s:addressLocality: "Cincinnati"
      s:addressRegion: "OH"
      s:postalCode: "45229"
      s:streetAddress: "3333 Burnet Ave"
      s:telephone: "+1(513)636-4200"
    s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
    s:department:
    - class: s:Organization
      s:legalName: "Allergy and Immunology"
      s:department:
      - class: s:Organization
        s:legalName: "Barski Research Lab"
        s:member:
        - class: s:Person
          s:name: Michael Kotliar
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/rgt-thor.md


doc: |
  What is THOR?
  --------------

  THOR is an HMM-based approach to detect and analyze differential peaks in two sets of ChIP-seq data
  from distinct biological conditions with replicates. THOR performs genomic signal processing, peak
  calling and p-value calculation in an integrated framework.

  For more information please refer to:
  -------------------------------------

  Allhoff, M., Sere K., Freitas, J., Zenke, M., Costa, I.G. (2016), Differential Peak Calling of ChIP-seq
  Signals with Replicates with THOR, Nucleic Acids Research, epub gkw680.