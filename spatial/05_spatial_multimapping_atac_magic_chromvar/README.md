# Spatial Multi-Mapping ATAC-MAGIC-ChromVAR Pipeline

This folder contains a three-step workflow for reconstructing spatial chromatin accessibility and inferring spatial regulatory programs.

## Overview

1. **`1.spATAC_reconstruction.ipynb`**  
   Reconstructs pseudo spatial-ATAC (`spATAC`) by integrating cross-modality mappings among snATAC, snRNA, and spatial transcriptomics.

2. **`2.spATAC_visualization.ipynb`**  
   Performs downstream visualization of reconstructed `spATAC`, including spatial embedding views, coordinate transformations, and locus-level plotting.

3. **`3.spATAC_chromvar.ipynb`**  
   Uses reconstructed `spATAC` to run chromVAR-based motif deviation analysis, estimate TF activity signals, and visualize spatial regulatory patterns.

## Data Flow

- **Input:** snATAC, snRNA, and spatial data with matched/transferable cell-state information.
- **Intermediate output:** reconstructed pseudo spatial-ATAC profiles (`spATAC`).
- **Final output:** spatially resolved TF activity/regulatory maps from chromVAR.

## Practical Notes

- Keep sample IDs and coordinate systems consistent across all three stages.
- Confirm sparse matrix handling and feature naming consistency before chromVAR scoring.
- If running large datasets, tune memory-related parameters (for example, chunk size, workers, and device selection).
