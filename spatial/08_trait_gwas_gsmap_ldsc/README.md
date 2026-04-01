# Spatial GWAS, Trait Interpretation, and LDSC Enrichment

This directory contains a five-step workflow for spatial trait association analysis in human spinal cord data, together with utility scripts and resource index files for reproducible GWAS processing.

## Script Overview

1. **`1.gsmap_gw17.sh`**  
   Runs GSMAP on the GW17 sample to generate spatial GWAS association outputs.

2. **`2.associations_analysis.ipynb`**  
   Performs downstream interpretation of spatial GWAS associations, including annotation-level summaries and trait-cell-type relationship analyses.

3. **`3.trait_heatmap.ipynb`**  
   Generates trait-level heatmaps from spatial GWAS outputs for cross-annotation visualization and comparison.

4. **`4.snATAC_ldscore.sh`**  
   Uses cell-type differential peaks derived from single-cell ATAC-seq to run LDSC-based trait enrichment analysis.

5. **`5.gwas_enrichment.ipynb`**  
   Visualizes enrichment results as heatmaps and related summary plots.

## Additional Utility Files

- **`format_sumstats_ALS_hg38.sh`**  
  Demonstrates how raw GWAS summary statistics are standardized into GSMAP-compatible input fields (e.g., `variant_id`, alleles, effect size, SE, P value, and sample size).

- **`gsmap_resource_structure.txt`**  
  Documents the directory structure of genome resource files used by the GSMAP workflow.

- **`ldsc_resource_structure.txt`**  
  Documents the directory structure of genome resource files used by the LDSC enrichment workflow.

## Workflow Logic

- **Step 1 (GSMAP):** compute spatially resolved GWAS association metrics for the GW17 sample.
- **Step 2-3 (Interpretation):** summarize and interpret the spatial GWAS signal across annotations, traits, and cell-type contexts.
- **Step 4 (LDSC):** quantify trait enrichment using snATAC differential CRE/peak-derived annotations.
- **Step 5 (Visualization):** present enrichment patterns in publication-ready heatmaps.

## Inputs and Outputs

- **Input:** spatial transcriptomics annotations, GWAS summary statistics, and snATAC differential peak annotations.
- **Intermediate outputs:** per-trait association tables and annotation-level significance matrices.
- **Final outputs:** trait enrichment summaries and heatmap visualizations.
