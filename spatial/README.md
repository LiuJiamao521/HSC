# Spatial Modules

This directory organizes Stereo-seq and cross-modal spatial analysis workflows.

- `01_stereo_preprocessing/`: GEM to analysis object conversion and basic input preparation.
- `02_integration_seurat_harmony_knn/`: stage-matched integration between single-cell and spatial data.
- `03_cell2location/`: spatial deconvolution and cell abundance estimation.
- `04_spatial_multimapping_atac_magic_chromvar/`: RNA-ATAC-spatial projection, MAGIC imputation, and chromVAR analysis.
- `05_nichecompass_niches/`: niche inference and niche-associated downstream analysis.
- `06_neighborhood_ncem_cellchat/`: neighborhood enrichment, NCEM, and CellChat communication analysis.
- `07_trait_gwas_gsmap_ldsc/`: trait enrichment and LDSC-related code collected from the current source tree.
