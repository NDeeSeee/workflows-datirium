cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement


'sd:upstream':
  chipseq_sample:
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-chipseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-atacseq-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  scores_files:
    type: File[]
    label: "ChIP/ATAC-Seq sample(s)"
    doc: "bigWig file from ChIP/ATAC-Seq sample(s)"
    'sd:upstreamSource': "chipseq_sample/bigwig"
    'sd:localLabel': true

  scores_labels:
    type: string[]
    label: "ChIP/ATAC-Seq sample(s)"
    doc: "Aliases for ChIP/ATAC-Seq sample(s)"
    'sd:upstreamSource': "chipseq_sample/alias"

  regions_files:
    type:
    - File
    - type: array
      items: File
    label: "BED file(s) with the peaks or genes"
    doc: |
      For a gene list that will be centered by
      TSS, create a file in a BED format with
      the columns that have the following
      meaning: "chrom start end name score strand".
      For a peak file that will be centered by
      the peak centers create a file in a BED
      format with the columns that have the
      following meaning: "chrom start end".

  before_region_start_length:
    type: int?
    default: 5000
    label: "Distance upstream of the reference-point selected"
    doc: |
      Distance upstream of the reference-point selected.
    'sd:layout':
      advanced: true
  
  after_region_start_length:
    type: int?
    default: 5000
    label: "Distance downstream of the reference-point selected"
    doc: |
      Distance downstream of the reference-point selected.
    'sd:layout':
      advanced: true

  bin_size:
    type: int?
    default: 10
    label: "Length, in bases, of the non-overlapping bins for averaging the score over the regions length"
    doc: |
      Length, in bases, of the non-overlapping bins for averaging the score over
      the regions length.
    'sd:layout':
      advanced: true

  plot_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "lines"
      - "fill"
      - "se"
      - "std"
    default: "lines"
    label: "Plot type to display"
    doc: |
      “lines” will plot the profile line based on the average type selected.
      “fill” fills the region between zero and the profile curve. The fill in
      color is semi transparent to distinguish different profiles. “se” and
      “std” color the region between the profile and the standard error or
      standard deviation of the data.
    'sd:layout':
      advanced: true

  sort_regions:
    type:
    - "null"
    - type: enum
      symbols:
      - "descend"
      - "ascend"
      - "no"
      - "keep"
    default: "descend"
    label: "Sorting order for regions"
    doc: |
      Whether the heatmap should present the regions sorted. The default is to sort in
      descending order based on the mean value per region. Note that “keep” and “no” are
      the same thing.
    'sd:layout':
      advanced: true

  what_to_show:
    type:
    - "null"
    - type: enum
      symbols:
      - plot, heatmap and colorbar
      - plot and heatmap
      - heatmap only
      - heatmap and colorbar
    default: "plot, heatmap and colorbar"
    label: "What show on the plot"
    doc: |
      The default is to include a summary or profile plot on top of the heatmap and a heatmap colorbar.
      Other options are: “plot and heatmap”, “heatmap only”, “heatmap and colorbar”, and the default “plot,
      heatmap and colorbar”.
    'sd:layout':
      advanced: true

  per_group:
    type: boolean?
    default: false
    label: "Plot all samples by group of regions instead of groups of regions by sample"
    doc: |
      The default is to plot all groups of regions by sample. Using this option instead plots all
      samples by group of regions. Note that this is only useful if you have multiple groups of
      regions by sample rather than group.
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 6
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  scores_matrix:
    type: File
    outputSource: compute_score_matrix/scores_matrix
    label: "Scores per genome regions matrix"
    doc: |
      Scores per genome regions matrix. This file that can be used
      with plotHeatmap and plotProfiles.

  heatmap_png:
    type: File
    outputSource: make_heatmap/heatmap_file
    label: "Heatmap for scores around centers of provided regions"
    doc: |
      Heatmap for scores around centers of provided regions.
      PNG format
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Heatmap for scores around centers of provided regions'

  compute_score_matrix_stdout_log:
    type: File
    outputSource: compute_score_matrix/stdout_log
    label: "compute_score_matrix stdout log"
    doc: |
      compute_score_matrix stdout log

  compute_score_matrix_stderr_log:
    type: File
    outputSource: compute_score_matrix/stderr_log
    label: "compute_score_matrix stderr log"
    doc: |
      compute_score_matrix stderr log

  make_heatmap_stdout_log:
    type: File
    outputSource: make_heatmap/stdout_log
    label: "make_heatmap stdout log"
    doc: |
      make_heatmap stdout log

  make_heatmap_stderr_log:
    type: File
    outputSource: make_heatmap/stderr_log
    label: "make_heatmap stderr log"
    doc: |
      make_heatmap stderr log


steps:

  recenter_regions:
    run:
      cwlVersion: v1.0
      class: Workflow
      requirements:
      - class: ScatterFeatureRequirement
      inputs:
        regions_files:
          type: File[]
      outputs:
        recentered_regions_files:
          type: File[]
          outputSource: recenter/output_file
      steps:
        recenter:
          run: ../tools/custom-bash.cwl
          in:
            input_file: regions_files
            script:
              default: |
                NUM_COL=$(head -n 1 "$0" | tr -d "\r" | tr "," "\t" | awk '{print NF; exit}')
                echo "Number of column in the regions file ${NUM_COL}"
                if [ "$NUM_COL" -ge 5 ]
                then
                  # BED for gene list
                  # chrom  start  end  name  [score] strand
                  echo "Recenter by the gene TSS"
                  cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{tss=$2; if ($6=="-") tss=$3; print $1"\t"tss"\t"tss+1}' > `basename $0`
                else
                  # BED for peaks
                  # chrom  start  end
                  echo "Recenter by the peak center"
                  cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{center=$2+int(($3-$2)/2); print $1"\t"center"\t"center+1}' > `basename $0`
                fi
          scatter: input_file
          out:
          - output_file
    in:
      regions_files: regions_files
    out:
    - recentered_regions_files

  compute_score_matrix:
    run: ../tools/deeptools-computematrix-referencepoint.cwl
    in:
      score_files: scores_files
      regions_files: recenter_regions/recentered_regions_files
      reference_point: 
        default: "TSS"  # doesn't matter what we set here because we centered regions ourlselves
      before_region_start_length: before_region_start_length
      after_region_start_length: after_region_start_length
      bin_size: bin_size
      sort_regions: sort_regions
      samples_label: scores_labels
      output_filename:
        default: "score_matrix.gz"
      missing_data_as_zero:
        default: true
      threads: threads
    out:
    - scores_matrix
    - stdout_log
    - stderr_log

  make_heatmap:
    run: ../tools/deeptools-plotheatmap.cwl
    in:
      plot_title:
        default: "Tag density around peak centers"
      scores_matrix: compute_score_matrix/scores_matrix
      output_filename:
        default: "score_heatmap.png"
      plot_type: plot_type
      sort_regions: sort_regions
      average_type_summary_plot:
        default: "mean"
      what_to_show: what_to_show
      ref_point_label:
        default: "Center/TSS"
      regions_label:
        source: recenter_regions/recentered_regions_files
        valueFrom: |
          ${
            var results = [];
            for (var i = 0; i < self.length; i++){
              results.push(self[i].basename.split('.').slice(0,-1).join('.'));
            }
            return results;
          }
      samples_label: scores_labels
      x_axis_label:
        default: "distance (bp)"
      y_axisLabel:
        default: "Signal mean"
      per_group: per_group
      plot_file_format:
        default: "png"
      legend_location:
        default: "upper-left"
    out:
    - heatmap_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "deepTools - heatmap for scores associated with genomic regions"
label: "deepTools - heatmap for scores associated with genomic regions"
s:alternateName: "Plots heatmap for scores around centers of provided regions"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/heatmap-deeptools.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  deepTools - heatmap for scores associated with genomic regions
  ======================================================

  Plots heatmap for scores around centers of provided regions