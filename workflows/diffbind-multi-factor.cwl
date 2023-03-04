cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_by_common_delim = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  dna_experiment:
  - "chipseq-se.cwl"
  - "chipseq-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/chipseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/chipseq-pe.cwl"
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
    label: "Experiment short name"
    sd:preview:
      position: 1

  alignment_files:
    type: File[]
    secondaryFiles:
    - .bai
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Sorted and indexed alignment files in bam format
    'sd:upstreamSource': "dna_experiment/bambai_pair"
    'sd:localLabel': true

  peak_files:
    type: File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc:
      Peak files in the MACS2 xls format. Number and order of the
      files should correspond to the files provided in --alignments
      parameter.
    'sd:upstreamSource': "dna_experiment/macs2_called_peaks"
    'sd:localLabel': true

  dataset_names:
    type: string[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Unique names for datasets provided in --alignments and --peaks
      parameters, no special characters or spaces are allowed. Number
      and order of the names should correspond to the values provided
      in --alignments and --peaks parameters.
    'sd:upstreamSource': "dna_experiment/alias"
    'sd:localLabel': true

  genome_coverage_files:
    type: File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Genome coverage files in bigWig format
    'sd:upstreamSource': "dna_experiment/bigwig"
    'sd:localLabel': true

  narrow_peak_files:
    type:
    - "null"
    - File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Called peaks files in narrowPeak format
    'sd:upstreamSource': "dna_experiment/macs2_narrow_peaks"
    'sd:localLabel': true

  broad_peak_files:
    type:
    - "null"
    - File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Called peaks files in broadPeak format
    'sd:upstreamSource': "dna_experiment/macs2_broad_peaks"
    'sd:localLabel': true

  annotation_file:
    type: File
    label: "Genome annotation file in TSV format"
    doc: |
      Genome annotation file in TSV format
    'sd:upstreamSource': "genome_indices/annotation"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    doc: |
      Chromosome length file
    'sd:upstreamSource': "genome_indices/chrom_length"

  design_formula:
    type: string
    label: "Design formula comprised of the metadata columns names"
    doc: |
      Design formula comprised of the metadata columns names.
      It should start with ~

  contrast:
    type: string?
    default: null
    label: "Contrast applied to the analysis results when calculating log2 fold changes"
    doc: |
      Contrast applied to the analysis results when calculating log2 fold changes.
      It should be formatted as a mathematical formula of values present in the
      metadata table. It is a required parameter if --method is set to edger. If not
      provided and --method is set to deseq2, the last term from the design formula
      will be used.

  base_levels:
    type: string?
    default: null
    label: "Base levels for each of the metadata columns"
    doc: |
      Base levels for each of the metadata columns. Number and order of the provided
      values should correspond to the metadata columns. Default: define base levels
      alphabetically.

  overlap_threshold:
    type: int?
    default: 2
    label: "Filtering threshold to keep only those peaks that are present in at least this many datasets"
    doc: |
      Filtering threshold to keep only those peaks that are present in at
      least this many datasets when generating consensus set of peaks.
      Ignored if --groupby is provided.
      Default: 2

  groupby:
    type: string?
    default: null
    label: "Column(s) from the metadata table to define datasets groups for obtaining the common peaks within each of them"
    doc: |
      Column(s) from the metadata table to define datasets groups for obtaining
      the common peaks within each of them. Union of such common peaks will be
      used as consensus peaks.
      Default: do not search for common peaks, use --minoverlap parameter instead

  metadata_file:
    type: File
    label: "TSV/CSV metadata file to describe datasets"
    doc: |
      TSV/CSV metadata file to describe datasets provided in --alignments
      and --peaks parameters. First column should have the name 'sample',
      all other columns names should be selected from the following list:
      Tissue, Factor, Condition, Treatment, Caller, Replicate. The values
      from the 'sample' column should correspond to the values provided in
      --aliases parameter. For a proper --contrast intepretation, values
      defined in each metadata column should not be used in any of the other
      columns. All metadata columns are treated as factors (no covariates
      are supported).

  scoreby:
    type:
    - "null"
    - type: enum
      symbols:
      - "pvalue"
      - "qvalue"
    default: "pvalue"
    label: "Score metrics to exclude low quality peaks"
    doc: |
      Score metrics to build peak overlap correlation heatmap and exclude low
      quality peaks based on the threshold provided in --score parameter.
      Default: pvalue
    'sd:layout':
      advanced: true

  score_threshold:
    type: float?
    default: 0.05
    label: "Filtering threshold for pvalue/qvalue"
    doc: |
      Filtering threshold to keep only those peaks where the metric selected
      in --scoreby parameter is less than or equal to the provided value.
      Default: 0.05
    'sd:layout':
      advanced: true

  rpkm_threshold:
    type: float?
    default: 1
    label: "Filtering threshold for maximum RPKM"
    doc: |
      Filtering threshold to keep only those peaks where the max RPKM for
      all datasets is bigger than or equal to the provided value.
      Default: 1
    'sd:layout':
      advanced: true

  padj_threshold:
    type: float?
    default: 0.05
    label: "Filtering threshold to report only significant differentially bound sites"
    doc: |
      Filtering threshold to report only differentially bound sites with adjusted
      P-value less than or equal to the provided value.
      Default: 0.05
    'sd:layout':
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: |
      Maximum distance from gene TSS (in both direction) overlapping which
      the peak will be assigned to the promoter region.
      Default: 1000 bp
    'sd:layout':
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: |
      Maximum distance from the promoter (only in upstream direction) overlapping
      which the peak will be assigned to the upstream region. Default: 20,000 bp"
    'sd:layout':
      advanced: true

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
      - "none"
    default: "none"
    label: "Hopach clustering method to be run on normalized read counts"
    doc: |
      Hopach clustering method to be run on normalized read counts.
      Default: do not run clustering
    'sd:layout':
      advanced: true

  row_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "cosangle"
    label: "Distance metric for HOPACH row clustering"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided.
      Default: cosangle
    'sd:layout':
      advanced: true

  column_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "euclid"
    label: "Distance metric for HOPACH column clustering"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided.
      Default: euclid
    'sd:layout':
      advanced: true

  center_row:
    type: boolean?
    default: false
    label: "Apply mean centering for normalized read counts prior to running clustering by row"
    doc: |
      Apply mean centering for normalized read counts prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
      Default: 1
    'sd:layout':
      advanced: true


outputs:

  gc_files:
    type: File[]
    label: "Genome coverage"
    doc: |
      Genome coverage files in bigWig format
    outputSource: pipe/gc_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Genome coverage"
        height: 120

  np_files:
    type:
    - "null"
    - File[]
    label: "Called peaks (narrowPeak format)"
    doc: |
      Called peaks files in narrowPeak format
    outputSource: pipe/np_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks"
        displayMode: "COLLAPSE"
        height: 40

  bp_files:
    type:
    - "null"
    - File[]
    label: "Called peaks (broadPeak format)"
    doc: |
      Called peaks files in broadPeak format
    outputSource: pipe/bp_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks"
        displayMode: "COLLAPSE"
        height: 40

  diff_sts_bigbed:
    type: File
    label: "Differentially bound sites (bigBed format)"
    doc: |
      Differentially bound sites in bigBed format
    outputSource: bed_to_bigbed/bigbed_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        format: 'bigbed'
        name: "Differentially bound sites"
        height: 40

  pk_venn_plot_png:
    type: File?
    label: "Consensus peaks venn diagram"
    doc: |
      Consensus peaks venn diagram
      PNG format
    outputSource: diffbind/pk_venn_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Consensus peaks venn diagram'

  pk_vrlp_plot_png:
    type: File?
    label: "Peakset overlap rate"
    doc: |
      Peakset overlap rate
      PNG format
    outputSource: diffbind/pk_vrlp_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Peakset overlap rate'

  pk_scr_corr_plot_png:
    type: File?
    label: "Datasets correlation (peak score)"
    doc: |
      Datasets correlation (peak score)
      PNG format
    outputSource: diffbind/pk_scr_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Datasets correlation (peak score)'

  rw_rds_corr_plot_png:
    type: File?
    label: "Datasets correlation (raw reads)"
    doc: |
      Datasets correlation (raw reads)
      PNG format
    outputSource: diffbind/rw_rds_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Datasets correlation (raw reads)'

  nr_rds_corr_plot_png:
    type: File?
    label: "Datasets correlation (normalized reads)"
    doc: |
      Datasets correlation (normalized reads)
      PNG format
    outputSource: diffbind/nr_rds_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Datasets correlation (normalized reads)'

  pk_prfl_plot_png:
    type: File?
    label: "Peak profiles"
    doc: |
      Peak profiles
      PNG format
    outputSource: diffbind/pk_prfl_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'Peak profiles'

  diff_vlcn_plot_png:
    type: File?
    label: "Volcano plot for differentially bound sites"
    doc: |
      Volcano plot for differentially bound sites
      PNG format
    outputSource: diffbind/diff_vlcn_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'Volcano plot for differentially bound sites'

  diff_ma_plot_png:
    type: File?
    label: "MA-plot for differentially bound sites"
    doc: |
      MA-plot for differentially bound sites
      PNG format
    outputSource: diffbind/diff_ma_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'MA-plot for differentially bound sites'

  nr_rds_pca_1_2_plot_png:
    type: File?
    label: "PCA (1,2) of not filtered normalized counts"
    doc: |
      PCA (1,2) of not filtered normalized counts
      PNG format
    outputSource: diffbind/nr_rds_pca_1_2_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'PCA (1,2) of not filtered normalized counts'

  nr_rds_pca_2_3_plot_png:
    type: File?
    label: "PCA (2,3) of not filtered normalized counts"
    doc: |
      PCA (2,3) of not filtered normalized counts
      PNG format
    outputSource: diffbind/nr_rds_pca_2_3_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'PCA (2,3) of not filtered normalized counts'

  nr_rds_mds_html:
    type: File?
    outputSource: diffbind/nr_rds_mds_html
    label: "MDS plot of normalized counts"
    doc: |
      MDS plot of normalized counts.
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  diff_sts_tsv:
    type: File
    label: "Differentially bound sites with assigned nearest genes"
    doc: |
      Differentially bound sites with assigned nearest genes
      TSV format
    outputSource: restore_columns/output_file
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Differentially bound sites'
          Title: 'Differentially bound sites'

  diff_sts_labeled_tsv:
    type: File
    label: "Differentially bound sites with labels"
    doc: |
      Differentially bound sites with labels
      TSV format
    outputSource: add_label_column/output_file

  volcano_plot_html_file:
    type: File
    label: "Volcano Plot"
    doc: |
      HTML index file for Volcano Plot
    outputSource: make_volcano_plot/html_file
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for Volcano Plot
    outputSource: make_volcano_plot/html_data

  ma_plot_html_file:
    type: File
    label: "MA-plot"
    doc: |
      HTML index file for MA-plot
    outputSource: make_ma_plot/html_file
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  ma_plot_html_data:
    type: Directory
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for MA-plot
    outputSource: make_ma_plot/html_data

  heatmap_html:
    type: File
    label: "Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format
    outputSource: morpheus_heatmap/heatmap_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  diffbind_stdout_log:
    type: File
    label: "DiffBind stdout log"
    doc: |
      DiffBind stdout log
    outputSource: diffbind/stdout_log

  diffbind_stderr_log:
    type: File
    label: "DiffBind stderr log"
    doc: |
      DiffBind stderr log
    outputSource: diffbind/stderr_log

  morpheus_stdout_log:
    type: File
    label: "Morpheus heatmap stdout log"
    doc: |
      Morpheus heatmap stdout log
    outputSource: morpheus_heatmap/stdout_log

  morpheus_stderr_log:
    type: File
    label: "Morpheus heatmap stderr log"
    doc: |
      Morpheus heatmap stderr log
    outputSource: morpheus_heatmap/stderr_log


steps:

  pipe:
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      inputs:
        genome_coverage_files:
          type: File[]
        narrow_peak_files:
          type:
          - "null"
          - File[]
        broad_peak_files:
          type:
          - "null"
          - File[]
      outputs:
        gc_files:
          type: File[]
        np_files:
          type:
          - "null"
          - File[]
        bp_files:
          type:
          - "null"
          - File[]
      expression: |
        ${
          var results = {};
          var output_names = [
            "gc_files",
            "np_files",
            "bp_files"
          ];
          var sources = [
            inputs.genome_coverage_files,
            inputs.narrow_peak_files,
            inputs.broad_peak_files
          ];
          for (var i = 0; i < sources.length; i++){
            var current_source = sources[i];
            var current_output_name = output_names[i];
            results[current_output_name] = null;
            if (current_source != null && current_source.length > 0){
              for (var j = 0; j < current_source.length; j++){
                    var new_item = current_source[j];
                    new_item["basename"] = "u" + "_" + i + "_" + j+ "_" + new_item.basename;
                    if (results[current_output_name] == null){
                      results[current_output_name] = [new_item];
                    } else {
                      results[current_output_name].push(new_item);
                    }
              }
            }
          }
          return results;
        }
    in:
      genome_coverage_files: genome_coverage_files
      narrow_peak_files: narrow_peak_files
      broad_peak_files: broad_peak_files
    out:
    - gc_files
    - np_files
    - bp_files

  diffbind:
    run: ../tools/diffbind-multi-factor.cwl
    in:
      alignment_files_: alignment_files
      peak_files: peak_files
      dataset_names: dataset_names
      metadata_file: metadata_file
      scoreby: scoreby
      score_threshold: score_threshold
      rpkm_threshold: rpkm_threshold
      overlap_threshold: overlap_threshold
      groupby:
        source: groupby
        valueFrom: $(split_by_common_delim(self))
      design_formula: design_formula
      contrast:
        source: contrast
        valueFrom: $(self==""?null:self)                 # safety measure
      base_levels:
        source: base_levels
        valueFrom: $(split_by_common_delim(self))
      analysis_method:
        default: "deseq2"                                # hardcoded to always use DESeq2 because EdgeR fails to run without contrast
      normalization_method:
        default: "auto"                                  # harcoded to auto as we don't allow to use EdgeR
      padj_threshold: padj_threshold
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row: center_row
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - pk_venn_plot_png
    - pk_vrlp_plot_png
    - pk_scr_corr_plot_png
    - rw_rds_corr_plot_png
    - nr_rds_corr_plot_png
    - pk_prfl_plot_png
    - diff_vlcn_plot_png
    - diff_ma_plot_png
    - nr_rds_pca_1_2_plot_png
    - nr_rds_pca_2_3_plot_png
    - nr_rds_mds_html
    - diff_sts_tsv
    - nr_rds_gct
    - stdout_log
    - stderr_log

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diffbind/diff_sts_tsv
      script:
        default: >
          cat $0 | grep -v "Start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  assign_genes:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: filter_columns/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
      - assign_genes/result_file
      - diffbind/diff_sts_tsv
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > iaintersect_result.tsv
          cat $1 | grep -v "Start" > diffbind_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "Refseq_id\tGene_id\ttxStart\ttxEnd\tStrand\tRegion\t${HEADER}" > `basename $0`;
          cat iaintersect_result.tsv | paste - diffbind_result.tsv >> `basename $0`
          rm iaintersect_result.tsv diffbind_result.tsv
    out:
    - output_file

  convert_to_bed:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      script:
        default: |
          cat "$0" | awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} } NR>1 {color="255,0,0"; if ($ix["log2FoldChange"]<0) color="0,255,0"; print $ix["Chr"]"\t"$ix["Start"]"\t"$ix["End"]"\tpvalue="$ix["pvalue"]";padj="$ix["padj"]";log2FC="$ix["log2FoldChange"]"\t"1000"\t"$ix["Strand"]"\t"$ix["Start"]"\t"$ix["End"]"\t"color}' > `basename $0`
    out:
    - output_file
  sort_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: convert_to_bed/output_file
      key:
        default: ["1,1","2,2n"]
    out:
    - sorted_file

  bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: sort_bed/sorted_file
      bed_type:
        default: "bed4+5"
      chrom_length_file: chrom_length_file
      output_filename:
        source: sort_bed/sorted_file
        valueFrom: $(self.basename.split('.').slice(0,-1).join('.') + ".bigBed")
    out: 
    - bigbed_file

  add_label_column:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diffbind/diff_sts_tsv
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > diff_sts_labeled.tsv;
          cat "$0" | awk -F "\t" '{print $1":"$2"-"$3"\t"$0}' >> diff_sts_labeled.tsv
    out:
    - output_file

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  make_ma_plot:
    run: ../tools/ma-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "baseMean"
      y_axis_column:
        default: "log2FoldChange"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: diffbind/nr_rds_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "DiffBind Multi-factor Analysis"
s:name: "DiffBind Multi-factor Analysis"
s:alternateName: "Runs DiffBind multi-factor analysis with manual control over major parameters"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/diffbind-multi-factor.cwl
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
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898


doc: |
  DiffBind Multi-factor Analysis
  ------------------------------
  
  DiffBind processes ChIP-Seq data enriched for genomic loci where specific protein/DNA binding occurs, including peak sets identified by ChIP-Seq peak callers and
  aligned sequence read datasets. It is designed to work with multiple peak sets simultaneously, representing different ChIP experiments (antibodies, transcription
  factor and/or histone marks, experimental conditions, replicates) as well as managing the results of multiple peak callers.

  For more information please refer to:
  -------------------------------------
  Ross-Innes CS, Stark R, Teschendorff AE, Holmes KA, Ali HR, Dunning MJ, Brown GD, Gojis O, Ellis IO, Green AR, Ali S, Chin S, Palmieri C, Caldas C, Carroll JS (2012).
  “Differential oestrogen receptor binding is associated with clinical outcome in breast cancer.” Nature, 481, -4.