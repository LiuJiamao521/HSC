# snRNA: Scanpy QC, annotation, and differential expression

**Scanpy** workflows that load Cell Ranger **gene expression** matrices into a multi-sample `AnnData`, run QC and embedding, integrate **metadata and cell-type labels**, and perform **differential expression** between groups or types.

## Workflow

1. **`1.Scanpy.ipynb`** (Part 1)  
   Read per-sample `filtered_feature_bc_matrix`, merge objects, and run standard Scanpy steps (QC, normalization, HVGs, neighbors, UMAP, clustering).

2. **`2.snRNA_annotation.ipynb`** (Part 2)  
   Import **metadata and cell-type labels** from CSV or related tables into `adata.obs`, aligned with the integrated object.

3. **`3.Gene_expression.ipynb`** (Part 3)  
   After annotation, run **differential expression** (e.g. `rank_genes_groups`), visualize markers, and export result tables.

## Input and output

- **Input:** Cell Ranger GEX `filtered_feature_bc_matrix` (optional annotation tables).  
- **Output:** Annotated `h5ad`, DE results, and figures; RNA-side input for velocity, GRN, CellOracle, etc.
