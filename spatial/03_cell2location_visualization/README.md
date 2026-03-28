# 03_cell2location_visualization

This directory contains scripts for **spot-level** spatial deconvolution and visualization, as well as a visualization script for single-cell mapping annotations generated upstream.

## Scripts

- `1.hsc_c2l.ipynb`  
  Core `cell2location` deconvolution workflow. It estimates cell-type abundance/composition for each spatial spot.

- `2.decon_visual.ipynb`  
  Visualization of the deconvolution results produced by `1.hsc_c2l.ipynb` (spot level), including spatial patterns of cell-type abundance.

- `x.mapping_visual.ipynb`  
  Visualization of single-cell mapping annotations generated from the `02_integration_seurat_harmony_knn` workflow.

## Workflow and Relationships

1. Run `1.hsc_c2l.ipynb` first to generate spot-level `cell2location` deconvolution results.  
2. Run `2.decon_visual.ipynb` to visualize the outputs from script 1.  
3. Use `x.mapping_visual.ipynb` independently to visualize mapping annotations from `02_integration_seurat_harmony_knn`.

## Notes

- Both `1.hsc_c2l.ipynb` and `2.decon_visual.ipynb` focus on **spot-level deconvolution**.  
- `x.mapping_visual.ipynb` focuses on **mapping-annotation visualization** and does not run `cell2location` deconvolution.
