cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type:
    - File
    - type: array
      items: File
    label: "RDS file(s) produced by BD Rhapsody Sequence Analysis Pipeline"
    doc: |
      Path to the RDS file(s) to load Seurat object(s)
      from. These files should be generated by the BD
      Rhapsody Sequence Analysis Pipeline and include
      Sample_Tag metadata column.

  sample_tags_metadata:
    type: File
    label: "Sample tags metadata file (first column should be sample_tag)"
    doc: |
      Path to the TSV/CSV file to assign names to the
      sample tags. This file must include exactly two
      columns. First column is a 'sample_tag'. It should
      correspond to all unique values from the 'Sample_Tag'
      column of the loaded Seurat object(s). Second column
      may have an arbitrary name. But it should include
      unique names for each sample tag.

  split_by_origin:
    type: boolean?
    default: false
    label: "Split each sample tag by origin"
    doc: |
      When assigning names to the sample tags, split 
      each of them by origin (the RDS file the data
      was loaded from). 
      Default: do not split
    "sd:layout":
      advanced: true

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
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
      - "5"
      - "6"
    default: "6"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 6
    "sd:layout":
      advanced: true


outputs:

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: sc_rna_load_rhapsody/feature_bc_matrices_folder
    label: "TAR-gzipped folder with the feature-barcode matrix in MEX format"
    doc: |
      Compressed folder with the merged feature-barcode 
      matrix from the loaded RDS files produced by the 
      BD Rhapsody Sequence Analysis Pipeline. 
      MEX format (TAR-gzipped).

  aggregation_metadata:
    type: File?
    outputSource: sc_rna_load_rhapsody/aggregation_metadata
    label: "Aggregation metadata in TSV format"
    doc: |
      Aggregation metadata file with names assigned
      to sample tags. The row order corresponds to
      the numeric suffixes added to cell barcodes in
      the merged feature-barcode matrix.
      TSV format.

  seurat_data_rds:
    type: File
    outputSource: sc_rna_load_rhapsody/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  sc_report_html_file:
    type: File?
    outputSource: sc_rna_load_rhapsody/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_rna_load_rhapsody_stdout_log:
    type: File
    outputSource: sc_rna_load_rhapsody/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_rna_load_rhapsody step.

  sc_rna_load_rhapsody_stderr_log:
    type: File
    outputSource: sc_rna_load_rhapsody/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_rna_load_rhapsody step.


steps:

  sc_rna_load_rhapsody:
    run: ../tools/sc-rna-load-rhapsody.cwl
    in:
      query_data_rds: query_data_rds
      sample_tags_metadata: sample_tags_metadata
      split_by_origin: split_by_origin
      export_html_report: export_html_report
      verbose:
        default: true
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - feature_bc_matrices_folder
    - aggregation_metadata
    - seurat_data_rds
    - sc_report_html_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell RNA-Seq BD Rhapsody Import"
s:name: "Single-cell RNA-Seq BD Rhapsody Import"
s:alternateName: "Single-cell RNA-Seq BD Rhapsody Import"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/sc-rna-load-rhapsody.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Single-cell RNA-Seq BD Rhapsody Import

  Imports RDS files produced by the BD Rhapsody
  Sequence Analysis Pipeline. Assigns names to the
  sample tags based on the provided metadata file.
  Exports results into compressed feature barcode
  matrix and RDS files.