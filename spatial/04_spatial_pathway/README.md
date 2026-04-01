# Spatial Pathway Scoring and Visualization

This directory contains a three-step workflow for pathway-level activity scoring and spatial visualization in spatial transcriptomics data.

## Workflow

1. **`1.ucell.ipynb`**  
   Computes pathway (gene-set) activity scores using UCell for predefined gene sets.

2. **`2.visual_pathway.ipynb`**  
   Visualizes pathway/UCell scores in spatial coordinates with publication-ready spatial maps.

3. **`3.spatial_expression.ipynb`**  
   Plots spatial expression patterns of selected genes across samples.

## Input and Output

- **Input:** spatial expression objects and predefined pathway gene sets.
- **Intermediate output:** per-cell or per-spot UCell score matrices for each pathway.
- **Final output:** spatial pathway-score figures and spatial gene-expression plots.

## Notes

- Keep gene symbols consistent with the expression matrix feature naming.
- Use consistent sample labels and coordinate systems across all three scripts.
- For large datasets, tune parallel settings and plotting parameters (for example, spot size and output resolution).
