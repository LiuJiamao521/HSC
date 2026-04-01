# Spatial Analysis Workflows

This directory contains the full spatial transcriptomics analysis framework used in this project, spanning preprocessing, cross-modal integration, pathway scoring, regulatory inference, niche modeling, intercellular communication, and trait-level GWAS enrichment.

## Module Overview

1. **`01_stereo_preprocessing/`**  
   Converts raw Stereo-seq `.tissue.gem.gz` files into analysis-ready Seurat objects, supports manual/optional niche-guided region selection, and exports filtered objects to `h5ad`.

2. **`02_integration_seurat_harmony_knn/`**  
   Performs integration between single-cell and spatial data, including annotation transfer, snRNA-snATAC linking, one-to-one snRNA-spatial mapping, tri-modal visualization (`snATAC -> snRNA -> spatial`), and conversion to `h5ad`.

3. **`03_cell2location_visualization/`**  
   Runs spot-level deconvolution with Cell2location and visualizes estimated cell-type abundance patterns in tissue space.

4. **`04_spatial_pathway/`**  
   Computes UCell-based pathway/gene-set scores, visualizes spatial pathway score distributions, and plots spatial gene-expression maps.

5. **`05_spatial_multimapping_atac_magic_chromvar/`**  
   Reconstructs pseudo spatial-ATAC (`spATAC`) from snATAC/snRNA/spatial mappings, performs downstream visualization, and infers spatial TF regulatory activity with chromVAR.

6. **`06_nichecompass_niches/`**  
   Identifies spatial niches from multi-sample spatial data and performs niche-level GO term enrichment for functional interpretation.

7. **`07_neighborhood_ncem_cellchat/`**  
   Characterizes local tissue architecture and communication via neighborhood enrichment (Squidpy), unsupervised NCEM, CellChat v2 with spatial constraints (multi-sample and single-sample), and optional non-spatial CellChat.

8. **`08_trait_gwas_gsmap_ldsc/`**  
   Performs trait association and enrichment analyses, including GSMAP on GW17 spatial data, downstream trait interpretation/heatmaps, snATAC differential-peak LDSC enrichment, enrichment visualization, GWAS sumstats formatting utilities, and GSMAP/LDSC resource index files.

## Recommended End-to-End Order

1. Preprocess Stereo-seq data (`01`).
2. Integrate with single-cell modalities and build cross-modal mappings (`02`).
3. Run spot deconvolution and visualization (`03`).
4. Compute and visualize pathway activity (`04`).
5. Build the spatial ATAC regulatory layer (`05`).
6. Infer niches and niche functions (`06`).
7. Model neighborhood-level communication (`07`).
8. Perform trait-level spatial GWAS and LDSC enrichment (`08`).

## Data Products Across the Spatial Pipeline

- **Core objects:** Seurat `rds` and AnnData `h5ad` spatial objects.
- **Cross-modal links:** snATAC-snRNA and snRNA-spatial mapping tables.
- **Pathway outputs:** UCell pathway score matrices and spatial pathway/gene-expression maps.
- **Regulatory outputs:** chromVAR TF activity maps and spatial regulatory patterns.
- **Niche outputs:** niche assignments and GO enrichment summaries.
- **Communication outputs:** NCEM interaction statistics and CellChat pathway networks.
- **Trait outputs:** GSMAP association matrices, LDSC enrichment tables, and publication-ready heatmaps.

## Notes

- Each module has its own `README.md` with script-level details and execution notes.
- Keep sample IDs, coordinate systems, and annotation names consistent across modules.
- For reproducibility, track software versions, thresholds, and resource paths used in each step.
