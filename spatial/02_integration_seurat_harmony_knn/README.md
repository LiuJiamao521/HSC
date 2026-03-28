# Spatial Integration, Annotation, and Tri-modal Mapping

This section describes two related but distinct analysis tasks:

1. **Spatial annotation transfer**: transferring cell-type annotations from single-cell RNA-seq to spatial transcriptomics data.
2. **Tri-modal mapping**: linking `snATAC`, `snRNA`, and spatial transcriptomics at the single-cell or single-spot level.

In this directory, **Script 1** completes the spatial annotation task, whereas **Scripts 2-4** together complete the tri-modal mapping task. **Script 5** converts the annotated spatial object into `h5ad` format for downstream Python-based analyses.

## Overview of the Workflow

1. Use `1.integration_rna_st.ipynb` to integrate single-cell RNA-seq and spatial transcriptomics data and transfer cell-type labels to spatial spots.
2. Use `2.link_rna_atac.ipynb` to link `snRNA` and `snATAC` profiles through shared cell barcodes from the same nucleus in the snMultiome data.
3. Use `3.mapping_rna_st.ipynb` to perform one-to-one mapping between `snRNA` cells and spatial cells/spots.
4. Use `4.draw_3d.ipynb` to visualize the links across `snATAC UMAP -> snRNA UMAP -> spatial physical coordinates`.
5. Use `5.convert.ipynb` to convert the annotated spatial `rds` object into `.h5ad`.

## Script 1: `1.integration_rna_st.ipynb`

### Purpose

This script transfers cell-type annotations from single-cell RNA-seq to spatial transcriptomics data. The output is a **spatial object with inferred cell-type annotations**.

### Key Idea

This step focuses on **annotation transfer**, not exact one-to-one cell pairing. In other words, the result is the predicted identity of each spatial spot or spatial cell, rather than a unique matched `snRNA` cell for every spatial observation.

### Main Functions

- Integrates single-cell RNA-seq and spatial transcriptomics data.
- Uses the single-cell reference annotation to infer cell identities in spatial data.
- Produces an annotated spatial object for downstream spatial analyses.

### Main Input

- Annotated single-cell RNA-seq reference.
- Preprocessed spatial transcriptomics object.

### Main Output

- Annotated spatial `rds` object.
- Spatial cell-type labels for each spatial observation.

## Script 2: `2.link_rna_atac.ipynb`

### Purpose

This script links `snRNA` and `snATAC` profiles using shared cell barcodes.

### Key Idea

Because the dataset is generated from **snMultiome**, RNA and ATAC are measured from the **same nucleus**. Therefore, this step is conceptually simple: the two modalities are linked directly through barcode identity.

### Main Functions

- Reads `snRNA` and `snATAC` metadata.
- Matches the two modalities by shared cell barcode.
- Generates a direct `snATAC <-> snRNA` correspondence table.

### Main Input

- `snRNA` metadata or embeddings.
- `snATAC` metadata or embeddings.
- Shared cell barcode information.

### Main Output

- `snATAC-snRNA` linkage table.
- Matched multiome cell pairs for downstream tri-modal analysis.

## Script 3: `3.mapping_rna_st.ipynb`

### Purpose

This script performs **one-to-one mapping** between `snRNA` cells and spatial cells/spots.

### Key Idea

This step is different from Script 1.

- **Script 1** transfers annotations from `snRNA` to spatial data, so the result is a **cell-type label** for each spatial spot.
- **Script 3** identifies a **paired relationship** between one `snRNA` cell and one spatial cell/spot, producing an explicit mapping table.

Therefore, Script 3 is used when a direct correspondence between individual `snRNA` and spatial observations is required for downstream multi-modal linking and visualization.

### Main Functions

- Integrates `snRNA` and spatial transcriptomics data in a shared space.
- Identifies matched `snRNA-spatial` pairs.
- Produces one-to-one mapping relationships for downstream tri-modal reconstruction.

### Main Input

- Single-cell RNA-seq object.
- Spatial transcriptomics object.
- Shared embedding or integrated representation.

### Main Output

- `snRNA-spatial` mapping table.
- Paired cell-level relationships between the two modalities.

## Script 4: `4.draw_3d.ipynb`

### Purpose

This script visualizes the tri-modal correspondence among `snATAC`, `snRNA`, and spatial transcriptomics.

### Key Idea

This step combines the two linkage results generated previously:

- `snATAC-snRNA` links from Script 2
- `snRNA-spatial` links from Script 3

By chaining these two relationships together, the script visualizes how a cell is represented across:

`snATAC UMAP -> snRNA UMAP -> spatial physical coordinates`

This provides an intuitive view of how chromatin accessibility, transcriptome state, and true tissue location are connected.

### Main Functions

- Reads `snATAC-snRNA` and `snRNA-spatial` linkage tables.
- Extracts embeddings from `snATAC` and `snRNA`.
- Extracts physical coordinates from spatial transcriptomics.
- Draws linked trajectories or correspondence lines across the three modalities.

### Main Input

- `snATAC-snRNA` links.
- `snRNA-spatial` links.
- `snATAC` UMAP coordinates.
- `snRNA` UMAP coordinates.
- Spatial coordinates.

### Main Output

- 3D or linked visualization of tri-modal mapping.
- Figures showing the connection from regulatory state to transcriptomic state to tissue position.

## Script 5: `5.convert.ipynb`

### Purpose

This script converts the annotated spatial `rds` object into `.h5ad` format.

### Main Functions

- Reads the annotated spatial Seurat object.
- Adds or preserves transferred cell-type annotations in metadata.
- Converts the object into AnnData format.
- Saves the output for downstream Python-based workflows.

### Main Input

- Annotated spatial `rds` object.

### Main Output

- Annotated spatial `.h5ad` object.

## Two Analysis Tasks in This Directory

### Task 1: Spatial annotation transfer

This task can be completed by **Script 1 alone**.

- Input: annotated single-cell RNA-seq reference and spatial transcriptomics data.
- Output: annotated spatial transcriptomics object.

### Task 2: Tri-modal mapping

This task requires **Scripts 2-4** together.

- Script 2 establishes the `snATAC-snRNA` relationship.
- Script 3 establishes the `snRNA-spatial` relationship.
- Script 4 visualizes the full `snATAC -> snRNA -> spatial` mapping chain.

## Recommended Execution Order

1. Run `1.integration_rna_st.ipynb` to generate spatial annotations.
2. Run `2.link_rna_atac.ipynb` to link `snATAC` and `snRNA` by barcode.
3. Run `3.mapping_rna_st.ipynb` to generate one-to-one `snRNA-spatial` mapping.
4. Run `4.draw_3d.ipynb` to visualize tri-modal mapping.
5. Run `5.convert.ipynb` to convert the annotated spatial object to `.h5ad`.

## Outputs of This Section

After completing this section, the following key outputs will be generated:

- An annotated spatial transcriptomics object.
- A direct `snATAC-snRNA` linkage table.
- A one-to-one `snRNA-spatial` mapping table.
- Visualization of `snATAC UMAP -> snRNA UMAP -> spatial coordinates`.
- An annotated spatial `.h5ad` object for downstream Python analysis.
