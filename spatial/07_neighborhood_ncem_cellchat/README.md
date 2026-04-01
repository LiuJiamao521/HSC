# Neighborhood, NCEM, and CellChat Analysis for Spatial Transcriptomics

This directory contains a three-part workflow to quantify local spatial architecture and intercellular communication in spinal cord spatial transcriptomics data.

## Scripts and Purpose

1. **`1.squidpy_ncem.ipynb`**  
   Performs neighborhood enrichment and unsupervised intercellular communication modeling with NCEM.

2. **CellChat V2 with spatial constraints**  
   - **`2.cellchat2_spatial_multisample.ipynb`**: joint analysis across multiple samples  
   - **`2.cellchat2_spatial_single_sample.ipynb`**: per-sample analysis  
   These notebooks quantify ligand-receptor communication while incorporating spatial distance constraints.

3. **`3.(option)cellchat_rna.ipynb`**  
   Performs conventional CellChat analysis without spatial constraints (treating data similarly to standard single-cell RNA workflows).

## Method Summary

### 1) Neighborhood enrichment and NCEM

To characterize spatial architecture and local cellular dependencies, neighborhood enrichment and communication modeling are performed within each sample:

- Neighborhood enrichment is computed with `squidpy.gr.nhood_enrichment` (Squidpy v1.6.1).
- Because the data are spot-level spatial transcriptomics, Cell2location (v0.1.4) is used to estimate cell-type abundance per spot.
- These spot-level cell-type compositions are used as input to NCEM (v0.1.4).
- A node-centric linear expression model is fit to predict gene expression from local cell-type composition and neighborhood context.
- Sender-receiver dependencies are constrained by a spatial connectivity graph derived from local spot neighborhoods.

Significant spatial interactions are retained using both criteria:

- Interaction magnitude (Euclidean norm of node-centric coefficients) `> 0.5`
- At least `50` differentially expressed genes with `q < 0.05` for sender-receiver-specific terms

Only significant communication links are visualized in circular network plots.

### 2) CellChat with spatial constraints

To identify specific molecular mediators of communication, CellChat (v2.1.2) is applied as a supervised ligand-receptor framework:

- CellChatDB.human is used as the curated ligand-receptor database.
- Communication probabilities are computed with `truncatedMean` (`trim = 0.1`).
- Spatial constraints are integrated:
  - `contact.range = 50` for contact-dependent signaling
  - `interaction.range = 100` for diffusible signaling
- Communication is summarized at the signaling-pathway level.

Filtering strategy:

- Keep interactions with `P < 0.05`
- Apply a top `30%` quantile filter on total interaction strength per pathway

Final weighted directed networks are visualized with chord plots, bubble plots, and heatmaps.

### 3) Conventional CellChat (no spatial constraint)

An optional baseline CellChat workflow is provided in `3.(option)cellchat_rna.ipynb`, where communication is inferred without explicit spatial-distance constraints.

## Input/Output

- **Input:** spatial transcriptomics expression matrix, spot coordinates, and cell-type composition (Cell2location-derived for NCEM workflow).
- **Output:** neighborhood enrichment statistics, NCEM-derived significant sender-receiver interactions, and CellChat pathway-level communication networks/figures.
