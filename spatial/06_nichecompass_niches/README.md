# NicheCompass Spatial Niche Analysis

This directory contains a two-step workflow for spatial niche discovery and functional interpretation.

## Workflow

1. **`1.stereo_seq_spinal_cord_all_sample.ipynb`**  
   Performs joint spatial clustering on all spatial transcriptomics samples to identify tissue niches and niche-associated spatial patterns.

2. **`2.niche_go.ipynb`**  
   Runs GO term enrichment analysis for genes associated with the identified niches and summarizes functional programs.

## Input and Output

- **Input:** integrated multi-sample spatial transcriptomics data with required metadata annotations.
- **Intermediate output:** niche assignments and niche-level gene sets.
- **Final output:** GO enrichment tables and visualization figures for niche-associated biological functions.
