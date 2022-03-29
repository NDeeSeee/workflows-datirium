cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: dummy_metadata.csv
    entry: |
      library_id
      Experiment
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $(inputs.vector_memory_limit * 1000000000)


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.1


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate
      experiment in MEX format. The rows consist of all the gene and peak features
      concatenated together and the columns are restricted to those barcodes that are
      identified as cells.

  aggregation_metadata:
    type: File?
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities.
      If --mex points to the Cell Ranger ARC Aggregate outputs, the aggr.csv
      file can be used. If Cell Ranger ARC Count outputs have been used in
      --mex, the file should include at least one column - 'library_id' and
      one row with the alias for Cell Ranger ARC Count experiment.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format. Tbi-index file is required.

  annotation_gtf_file:
    type: File
    inputBinding:
      prefix: "--annotations"
    doc: |
      Path to the genome annotation file in GTF format

  grouping_data:
    type: File?
    inputBinding:
      prefix: "--grouping"
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      Default: each dataset is assigned to a separate group.

  blacklisted_regions_file:
    type: File?
    inputBinding:
      prefix: "--blacklisted"
    doc: |
      Path to the blacklisted regions file in BED format

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells.
      Default: use all cells.

  gex_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--gexmincells"
    doc: |
      Include only GEX features detected in at least this many cells.
      Default: 5

  gex_minimum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--mingenes"
    doc: |
      Include cells where at least this many GEX features are detected.
      If multiple values provided, each of them will be applied to the
      correspondent dataset from the --mex input based on the --identity
      file.
      Default: 250 (applied to all datasets)

  gex_maximum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--maxgenes"
    doc: |
      Include cells with the number of GEX features not bigger than this value.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 5000 (applied to all datasets)

  gex_minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--gexminumi"
    doc: |
      Include cells where at least this many GEX UMIs (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 500 (applied to all datasets)

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial GEX features.
      Default: '^Mt-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of GEX transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5

  regress_mito_perc:
    type: boolean?
    inputBinding:
      prefix: "--regressmt"
    doc: |
      Regress mitochondrial genes expression as a confounding source of variation
      when identifying GEX based clusters for calling custom MACS2 peaks.
      Ignored if --callpeaks is not provided.
      Default: false

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      GEX as log10(genes)/log10(UMIs). If multiple values provided, each of them will
      be applied to the correspondent dataset from the --mex input based on the
      --identity file.
      Default: 0.8 (applied to all datasets)

  atac_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--atacmincells"
    doc: |
      Include only ATAC features detected in at least this many cells.
      Default: 5

  atac_minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--atacminumi"
    doc: |
      Include cells where at least this many ATAC UMIs (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 1000 (applied to all datasets)

  maximum_nucl_signal:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--maxnuclsignal"
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments. If multiple values provided, each of
      them will be applied to the correspondent dataset from the --mex input
      based on the --identity file
      Default: 4 (applied to all datasets)

  minimum_tss_enrich:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--mintssenrich"
    doc: |
      Include cells with the TSS enrichment score not lower than this value.
      Score is calculated based on the ratio of fragments centered at the TSS
      to fragments in TSS-flanking regions. If multiple values provided, each
      of them will be applied to the correspondent dataset from the --mex input
      based on the --identity file.
      Default: 2 (applied to all datasets)

  minimum_frip:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minfrip"
    doc: |
      Include cells with the FRiP not lower than this value. If multiple values
      provided, each of them will be applied to the correspondent dataset from
      the --mex input based on the --identity file.
      Default: 0.15 (applied to all datasets)

  maximum_blacklisted_ratio:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--maxblacklisted"
    doc: |
      Include cells with the ratio of fragments in genomic blacklist regions
      not bigger than this value. If multiple values provided, each of them
      will be applied to the correspondent dataset from the --mex input based
      on the --identity file.
      Default: 0.05 (applied to all datasets)

  call_peaks:
    type:
    - "null"
    - type: enum
      symbols:
      - "identity"
      - "cluster"
    inputBinding:
      prefix: "--callpeaks"
    doc: |
      Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.
      Peaks are called per identity (identity) or per GEX cluster (cluster) after applying
      all GEX related thresholds, maximum nucleosome signal, and minimum TSS enrichment
      score filters. If set to 'cluster' GEX clusters are identified based on the parameters
      set with --resolution, --gexndim, --highvargex, --gexnorm, and --skipgexntrg.
      Default: do not call peaks

  gex_normalization:
    type:
    - "null"
    - type: enum
      symbols:
      - "sct"
      - "log"
      - "sctglm"
    inputBinding:
      prefix: "--gexnorm"
    doc: |
      Normalization method to be used when identifying GEX based clusters for
      calling custom MACS2 peaks. Ignored if --callpeaks is not set to 'cluster'.
      Default: sct

  gex_high_var_features_count:
    type: int?
    inputBinding:
      prefix: "--highvargex"
    doc: |
      Number of highly variable GEX features to detect. Used for GEX datasets
      integration, scaling, and dimensional reduction when identifying GEX based
      clusters for calling custom MACS2 peaks. Ignored if --callpeaks is not set
      to 'cluster'.
      Default: 3000

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "seurat"
      - "none"
    inputBinding:
      prefix: "--ntgr"
    doc: |
      Integration method for GEX datasets when identifying GEX based clusters
      for calling custom MACS2 peaks. Automatically set to 'none' if --mex points
      to the Cell Ranger ARC Count outputs (single, not aggregated dataset that
      doesn't require any integration). Ignored if --callpeaks is not set to
      'cluster'.
      Default: seurat

  gex_dimensionality:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--gexndim"
    doc: |
      Dimensionality to use in GEX UMAP projection and clustering when identifying
      GEX based clusters for calling custom MACS2 peaks (from 1 to 50). If single
      number N is provided, use from 1 to N PCs. If multiple numbers are provided,
      subset to only selected PCs. Ignored if --callpeaks is not set to 'cluster'.
      Default: from 1 to 10

  resolution:
    type: float?
    inputBinding:
      prefix: "--resolution"
    doc: |
      Resolution to be used when identifying GEX based clusters for calling
      custom MACS2 peaks. Ignored if --callpeaks is not set to 'cluster'.
      Default: 0.3

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: "--h5seurat"
    doc: |
      Save Seurat data to h5seurat file.
      Default: false

  low_memory:
    type: boolean?
    inputBinding:
      prefix: "--lowmem"
    doc: |
      Attempts to minimize RAM usage when integrating multiple datasets
      with SCTransform algorithm (slows down the computation).
      Ignored if --callpeaks is not set to 'cluster', if --ntgr is not set
      to 'seurat', if --gexnorm is not set to either 'sct' or 'glm'.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./seurat

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  raw_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.png"
    doc: |
      Peak density per cell (not filtered).
      PNG format

  raw_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.pdf"
    doc: |
      Peak density per cell (not filtered).
      PDF format

  raw_bl_cnts_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_bl_cnts_dnst.png"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (not filtered).
      PNG format

  raw_bl_cnts_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_bl_cnts_dnst.pdf"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (not filtered).
      PDF format

  raw_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_pca_1_2_qc_mtrcs.png"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format

  raw_pca_1_2_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_pca_1_2_qc_mtrcs.pdf"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered).
      PDF format

  raw_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_pca_2_3_qc_mtrcs.png"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format

  raw_pca_2_3_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_pca_2_3_qc_mtrcs.pdf"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PDF format

  raw_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.png"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.pdf"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format

  raw_gex_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst.png"
    doc: |
      GEX UMI density per cell (not filtered).
      PNG format

  raw_gex_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst.pdf"
    doc: |
      GEX UMI density per cell (not filtered).
      PDF format

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.png"
    doc: |
      ATAC UMI density per cell (not filtered).
      PNG format

  raw_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.pdf"
    doc: |
      ATAC UMI density per cell (not filtered).
      PDF format

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Gene density per cell (not filtered).
      PNG format

  raw_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.pdf"
    doc: |
      Gene density per cell (not filtered).
      PDF format

  raw_gex_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gex_atac_umi_corr.png"
    doc: |
      GEX vs ATAC UMIs per cell correlation (not filtered).
      PNG format

  raw_gex_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gex_atac_umi_corr.pdf"
    doc: |
      GEX vs ATAC UMIs per cell correlation (not filtered).
      PDF format

  raw_tss_enrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_tss_enrch.png"
    doc: |
      TSS Enrichment Score (not filtered).
      PNG format

  raw_tss_enrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_tss_enrch.pdf"
    doc: |
      TSS Enrichment Score (not filtered).
      PDF format

  raw_frg_len_hist_png:
    type: File?
    outputBinding:
      glob: "*_raw_frg_len_hist.png"
    doc: |
      Fragments Length Histogram (not filtered).
      PNG format

  raw_frg_len_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_frg_len_hist.pdf"
    doc: |
      Fragments Length Histogram (not filtered).
      PDF format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PNG format

  raw_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PDF format

  raw_mito_perc_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst.png"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format

  raw_mito_perc_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst.pdf"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PDF format

  raw_nvlt_score_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst.png"
    doc: |
      Novelty score density per cell (not filtered).
      PNG format

  raw_nvlt_score_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst.pdf"
    doc: |
      Novelty score density per cell (not filtered).
      PDF format

  raw_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format

  raw_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (not filtered).
      PDF format


  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst.png"
    doc: |
      Peak density per cell (intermediate filtered).
      PNG format

  mid_fltr_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst.pdf"
    doc: |
      Peak density per cell (intermediate filtered).
      PDF format

  mid_fltr_bl_cnts_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_bl_cnts_dnst.png"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (intermediate filtered).
      PNG format

  mid_fltr_bl_cnts_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_bl_cnts_dnst.pdf"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (intermediate filtered).
      PDF format

  mid_fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_pca_1_2_qc_mtrcs.png"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PNG format

  mid_fltr_pca_1_2_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_pca_1_2_qc_mtrcs.pdf"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PDF format

  mid_fltr_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_pca_2_3_qc_mtrcs.png"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PNG format

  mid_fltr_pca_2_3_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_pca_2_3_qc_mtrcs.pdf"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PDF format

  mid_fltr_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_cell_count.png"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PNG format

  mid_fltr_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_cell_count.pdf"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PDF format

  mid_fltr_gex_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gex_umi_dnst.png"
    doc: |
      GEX UMI density per cell (intermediate filtered).
      PNG format

  mid_fltr_gex_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gex_umi_dnst.pdf"
    doc: |
      GEX UMI density per cell (intermediate filtered).
      PDF format

  mid_fltr_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst.png"
    doc: |
      ATAC UMI density per cell (intermediate filtered).
      PNG format

  mid_fltr_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst.pdf"
    doc: |
      ATAC UMI density per cell (intermediate filtered).
      PDF format

  mid_fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst.png"
    doc: |
      Gene density per cell (intermediate filtered).
      PNG format

  mid_fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst.pdf"
    doc: |
      Gene density per cell (intermediate filtered).
      PDF format

  mid_fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gex_atac_umi_corr.png"
    doc: |
      GEX vs ATAC UMIs per cell correlation (intermediate filtered).
      PNG format

  mid_fltr_gex_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gex_atac_umi_corr.pdf"
    doc: |
      GEX vs ATAC UMIs per cell correlation (intermediate filtered).
      PDF format

  mid_fltr_tss_enrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_enrch.png"
    doc: |
      TSS Enrichment Score (intermediate filtered).
      PNG format

  mid_fltr_tss_enrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_enrch.pdf"
    doc: |
      TSS Enrichment Score (intermediate filtered).
      PDF format

  mid_fltr_frg_len_hist_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frg_len_hist.png"
    doc: |
      Fragments Length Histogram (intermediate filtered).
      PNG format

  mid_fltr_frg_len_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frg_len_hist.pdf"
    doc: |
      Fragments Length Histogram (intermediate filtered).
      PDF format

  mid_fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (intermediate filtered).
      PNG format

  mid_fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (intermediate filtered).
      PDF format

  mid_fltr_mito_perc_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_perc_dnst.png"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered).
      PNG format

  mid_fltr_mito_perc_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_perc_dnst.pdf"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered).
      PDF format

  mid_fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_score_dnst.png"
    doc: |
      Novelty score density per cell (intermediate filtered).
      PNG format

  mid_fltr_nvlt_score_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_score_dnst.pdf"
    doc: |
      Novelty score density per cell (intermediate filtered).
      PDF format

  mid_fltr_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (intermediate filtered).
      PNG format

  mid_fltr_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (intermediate filtered).
      PDF format

  mid_ntgr_gex_elbow_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_ntgr_gex_elbow.png"
    doc: |
      Elbow plot from GEX PCA of filtered integrated/scaled datasets.
      PNG format

  mid_ntgr_gex_elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_ntgr_elbow.pdf"
    doc: |
      Elbow plot from GEX PCA of filtered integrated/scaled datasets.
      PDF format

  mid_ntgr_gex_qc_dim_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_ntgr_gex_qc_dim_corr.png"
    doc: |
      Correlation plots between main QC metrics and PCA reduction on GEX assay.
      PNG format

  mid_ntgr_gex_qc_dim_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_ntgr_gex_qc_dim_corr.pdf"
    doc: |
      Correlation plots between main QC metrics and PCA reduction on GEX assay.
      PDF format

  mid_clst_gex_umap_res_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_res_*.png"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PNG format

  mid_clst_gex_umap_res_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PDF format

  mid_clst_gex_umap_spl_by_cond_res_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_spl_by_cond_res_*.png"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered GEX datasets.
      PNG format

  mid_clst_gex_umap_spl_by_cond_res_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_spl_by_cond_res_*.pdf"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered GEX datasets.
      PDF format

  mid_clst_gex_umap_qc_mtrcs_res_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_qc_mtrcs_res_*.png"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered GEX datasets.
      PNG format

  mid_clst_gex_umap_qc_mtrcs_res_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_clst_gex_umap_qc_mtrcs_res_*.pdf"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered GEX datasets.
      PDF format


  fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst.png"
    doc: |
      Peak density per cell (filtered).
      PNG format

  fltr_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst.pdf"
    doc: |
      Peak density per cell (filtered).
      PDF format

  fltr_bl_cnts_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_bl_cnts_dnst.png"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (filtered).
      PNG format

  fltr_bl_cnts_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_bl_cnts_dnst.pdf"
    doc: |
      Density of fraction of fragments within blacklisted regions per cell (filtered).
      PDF format

  fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_pca_1_2_qc_mtrcs.png"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered).
      PNG format

  fltr_pca_1_2_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_pca_1_2_qc_mtrcs.pdf"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered).
      PDF format

  fltr_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_pca_2_3_qc_mtrcs.png"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PNG format

  fltr_pca_2_3_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_pca_2_3_qc_mtrcs.pdf"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PDF format

  fltr_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_cell_count.png"
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_cell_count.pdf"
    doc: |
      Number of cells per dataset (filtered).
      PDF format

  fltr_gex_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gex_umi_dnst.png"
    doc: |
      GEX UMI density per cell (filtered).
      PNG format

  fltr_gex_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gex_umi_dnst.pdf"
    doc: |
      GEX UMI density per cell (filtered).
      PDF format

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst.png"
    doc: |
      ATAC UMI density per cell (filtered).
      PNG format

  fltr_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst.pdf"
    doc: |
      ATAC UMI density per cell (filtered).
      PDF format

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst.png"
    doc: |
      Gene density per cell (filtered).
      PNG format

  fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst.pdf"
    doc: |
      Gene density per cell (filtered).
      PDF format

  fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gex_atac_umi_corr.png"
    doc: |
      GEX vs ATAC UMIs per cell correlation (filtered).
      PNG format

  fltr_gex_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gex_atac_umi_corr.pdf"
    doc: |
      GEX vs ATAC UMIs per cell correlation (filtered).
      PDF format

  fltr_tss_enrch_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_enrch.png"
    doc: |
      TSS Enrichment Score (filtered).
      PNG format

  fltr_tss_enrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_enrch.pdf"
    doc: |
      TSS Enrichment Score (filtered).
      PDF format

  fltr_frg_len_hist_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frg_len_hist.png"
    doc: |
      Fragments Length Histogram (filtered).
      PNG format

  fltr_frg_len_hist_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frg_len_hist.pdf"
    doc: |
      Fragments Length Histogram (filtered).
      PDF format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PNG format

  fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PDF format

  fltr_mito_perc_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_perc_dnst.png"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format

  fltr_mito_perc_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_perc_dnst.pdf"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (filtered).
      PDF format

  fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_score_dnst.png"
    doc: |
      Novelty score density per cell (filtered).
      PNG format

  fltr_nvlt_score_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_score_dnst.pdf"
    doc: |
      Novelty score density per cell (filtered).
      PDF format

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format

  fltr_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (filtered).
      PDF format


  seurat_filtered_data_rds:
    type: File
    outputBinding:
      glob: "*_filtered_data.rds"
    doc: |
      Filtered Seurat data in RDS format

  seurat_filtered_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_filtered_data.h5seurat"
    doc: |
      Filtered Seurat data in h5seurat format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_multiome_filter.R"]
arguments:
- valueFrom: |
    ${
      if (inputs.aggregation_metadata) {
        return inputs.aggregation_metadata;
      } else {
        return runtime.outdir + "/dummy_metadata.csv"
      }
    }
  prefix: "--identity"


stdout: sc_multiome_filter_stdout.log
stderr: sc_multiome_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Multiome Filter"
s:name: "Single-cell Multiome Filter"
s:alternateName: "Filters single-cell multiome datasets based on the common QC metrics"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-multiome-filter.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
  Single-cell Multiome Filter
  ===================
  Filters single-cell multiome datasets based on the common QC metrics


s:about: |
  usage: sc_multiome_filter.R
        [-h] --mex MEX --identity IDENTITY --fragments FRAGMENTS --annotations
        ANNOTATIONS [--grouping GROUPING] [--blacklisted BLACKLISTED]
        [--barcodes BARCODES] [--gexmincells GEXMINCELLS]
        [--mingenes [MINGENES ...]] [--maxgenes [MAXGENES ...]]
        [--gexminumi [GEXMINUMI ...]] [--mitopattern MITOPATTERN]
        [--maxmt MAXMT] [--regressmt] [--minnovelty [MINNOVELTY ...]]
        [--atacmincells ATACMINCELLS] [--atacminumi [ATACMINUMI ...]]
        [--maxnuclsignal [MAXNUCLSIGNAL ...]]
        [--mintssenrich [MINTSSENRICH ...]] [--minfrip [MINFRIP ...]]
        [--maxblacklisted [MAXBLACKLISTED ...]]
        [--callpeaks {identity,cluster}] [--gexnorm {sct,log,sctglm}]
        [--highvargex HIGHVARGEX] [--ntgr {seurat,none}]
        [--gexndim [GEXNDIM ...]] [--resolution RESOLUTION] [--pdf] [--verbose]
        [--h5seurat] [--lowmem] [--output OUTPUT] [--cpus CPUS]
        [--memory MEMORY]

  Seurat Multiome Filtering Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX             Path to the folder with feature-barcode matrix from
                          Cell Ranger ARC Count/Aggregate experiment in MEX
                          format. The rows consist of all the gene and peak
                          features concatenated together and the columns are
                          restricted to those barcodes that are identified as
                          cells.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities. If --mex points to the Cell Ranger ARC
                          Aggregate outputs, the aggr.csv file can be used. If
                          Cell Ranger ARC Count outputs have been used in --mex,
                          the file should include at least one column -
                          'library_id' and one row with the alias for Cell
                          Ranger ARC Count experiment.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          observed in the experiment in TSV format. Tbi-index
                          file is required.
    --annotations ANNOTATIONS
                          Path to the genome annotation file in GTF format
    --grouping GROUPING   Path to the TSV/CSV file to define datasets grouping.
                          First column - 'library_id' with the values provided
                          in the same order as in the correspondent column of
                          the --identity file, second column 'condition'.
                          Default: each dataset is assigned to a separate group.
    --blacklisted BLACKLISTED
                          Path to the optional blacklisted regions file in BED
                          format
    --barcodes BARCODES   Path to the headerless TSV/CSV file with the list of
                          barcodes to select cells of interest (one barcode per
                          line). Prefilters input feature-barcode matrix to
                          include only selected cells. Default: use all cells.
    --gexmincells GEXMINCELLS
                          Include only GEX features detected in at least this
                          many cells. Default: 5 (applied to all datasets)
    --mingenes [MINGENES ...]
                          Include cells where at least this many GEX features
                          are detected. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the --mex input based on the --identity file. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES ...]
                          Include cells with the number of GEX features not
                          bigger than this value. If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the --mex input based on the --identity
                          file. Default: 5000 (applied to all datasets)
    --gexminumi [GEXMINUMI ...]
                          Include cells where at least this many GEX UMIs
                          (transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the --mex input based on
                          the --identity file. Default: 500 (applied to all
                          datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial GEX features.
                          Default: '^Mt-'
    --maxmt MAXMT         Include cells with the percentage of GEX transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5 (applied to all datasets)
    --regressmt           Regress mitochondrial genes expression as a
                          confounding source of variation when identifying GEX
                          based clusters for calling custom MACS2 peaks. Ignored
                          if --callpeaks is not provided. Default: false
    --minnovelty [MINNOVELTY ...]
                          Include cells with the novelty score not lower than
                          this value, calculated for GEX as
                          log10(genes)/log10(UMIs). If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the --mex input based on the --identity
                          file. Default: 0.8 (applied to all datasets)
    --atacmincells ATACMINCELLS
                          Include only ATAC features detected in at least this
                          many cells. Default: 5 (applied to all datasets)
    --atacminumi [ATACMINUMI ...]
                          Include cells where at least this many ATAC UMIs
                          (transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the --mex input based on
                          the --identity file. Default: 1000 (applied to all
                          datasets)
    --maxnuclsignal [MAXNUCLSIGNAL ...]
                          Include cells with the nucleosome signal not bigger
                          than this value. Nucleosome signal quantifies the
                          approximate ratio of mononucleosomal to nucleosome-
                          free fragments. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the --mex input based on the --identity file Default:
                          4 (applied to all datasets)
    --mintssenrich [MINTSSENRICH ...]
                          Include cells with the TSS enrichment score not lower
                          than this value. Score is calculated based on the
                          ratio of fragments centered at the TSS to fragments in
                          TSS-flanking regions. If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the --mex input based on the --identity
                          file. Default: 2 (applied to all datasets)
    --minfrip [MINFRIP ...]
                          Include cells with the FRiP not lower than this value.
                          If multiple values provided, each of them will be
                          applied to the correspondent dataset from the --mex
                          input based on the --identity file. Default: 0.15
                          (applied to all datasets)
    --maxblacklisted [MAXBLACKLISTED ...]
                          Include cells with the ratio of fragments in genomic
                          blacklist regions not bigger than this value. If
                          multiple values provided, each of them will be applied
                          to the correspondent dataset from the --mex input
                          based on the --identity file. Default: 0.05 (applied
                          to all datasets)
    --callpeaks {identity,cluster}
                          Call peaks with MACS2 instead of those that are
                          provided by Cell Ranger ARC Count. Peaks are called
                          per identity (identity) or per GEX cluster (cluster)
                          after applying all GEX related thresholds, maximum
                          nucleosome signal, and minimum TSS enrichment score
                          filters. If set to 'cluster' GEX clusters are
                          identified based on the parameters set with
                          --resolution, --gexndim, --highvargex, --gexnorm, and
                          --skipgexntrg. Default: do not call peaks
    --gexnorm {sct,log,sctglm}
                          Normalization method to be used when identifying GEX
                          based clusters for calling custom MACS2 peaks. Ignored
                          if --callpeaks is not set to 'cluster'. Default: sct
    --highvargex HIGHVARGEX
                          Number of highly variable GEX features to detect. Used
                          for GEX datasets integration, scaling, and dimensional
                          reduction when identifying GEX based clusters for
                          calling custom MACS2 peaks. Ignored if --callpeaks is
                          not set to 'cluster'. Default: 3000
    --ntgr {seurat,none}  Integration method for GEX datasets when identifying
                          GEX based clusters for calling custom MACS2 peaks.
                          Automatically set to 'none' if --mex points to the
                          Cell Ranger ARC Count outputs (single, not aggregated
                          dataset that doesn't require any integration). Ignored
                          if --callpeaks is not set to 'cluster'. Default:
                          seurat
    --gexndim [GEXNDIM ...]
                          Dimensionality to use in GEX UMAP projection and
                          clustering when identifying GEX based clusters for
                          calling custom MACS2 peaks (from 1 to 50). If single
                          number N is provided, use from 1 to N PCs. If multiple
                          numbers are provided, subset to only selected PCs.
                          Ignored if --callpeaks is not set to 'cluster'.
                          Default: from 1 to 10
    --resolution RESOLUTION
                          Resolution to be used when identifying GEX based
                          clusters for calling custom MACS2 peaks. Ignored if
                          --callpeaks is not set to 'cluster'. Default: 0.3
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --lowmem              Attempts to minimize RAM usage when integrating
                          multiple datasets with SCTransform algorithm (slows
                          down the computation). Ignored if --callpeaks is not
                          set to 'cluster', if --ntgr is not set to 'seurat', if
                          --gexnorm is not set to either 'sct' or 'glm'.
                          Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32