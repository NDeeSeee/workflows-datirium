cwlVersion: v1.0
class: Workflow


inputs:

  genome_label:
    type: string?
    label: "Genome label"
    doc: "Genome label is used by web-ui to show label"
    'sd:preview':
      position: 1

  genome_description:
    type: string?
    label: "Genome description"
    doc: "Genome description is used by web-ui to show description"
    'sd:preview':
      position: 2

  genome_details:
    type: string?
    label: "Genome details"
    doc: "Genome details"
    'sd:preview':
      position: 3

  genome_folder:
    type: Directory
    label: "Genome folder"
    doc: "Genome folder with FASTA files"


outputs:

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Bismark generated indices folder"
    outputSource: prepare_indices/indices_folder

  genome_directory:
    type: Directory
    label: "Genome folder"
    doc: "Genome folder with FASTA files"
    outputSource: genome_folder


steps:

  prepare_indices:
    run: ../tools/bismark-prepare-genome.cwl
    in:
      genome_folder: genome_folder
    out: [indices_folder]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "bismark-indices"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/bismark-indices.cwl
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

doc: |
  Prepares indices for Bismark Methylation Analysis. Bowtie2 aligner is used by default.

s:about: |
  Runs bismark_genome_preparation script to prepare indices for Bismark Methylation Analysis. Bowtie2 aligner is used by default.
