# Spatial Preprocessing

This section describes how raw Stereo-seq data are processed from `.tissue.gem.gz` files into analysis-ready objects for downstream spatial transcriptomics workflows.

## Input Data

The raw input files for spatial transcriptomics are Stereo-seq output files in `.tissue.gem.gz` format. Each sample should first be organized with a consistent sample name before object conversion and region selection.

## Overall Workflow

1. Use Code 1 to convert raw `.tissue.gem.gz` files into Seurat `rds` objects.
2. Optionally use NicheCompass to perform unsupervised spatial clustering and identify biologically meaningful regions.
3. Use Code 2 to manually select target regions or spots in spatial coordinates, or directly select regions based on the unsupervised clustering results.
4. Convert the filtered Seurat `rds` objects to `h5ad` format for downstream Python/Scanpy analysis.

## Code 1: `1.gem2rds.ipynb` + `stomics_seurat_only_gem2rds.R`

This step converts raw `.tissue.gem.gz` files into Seurat spatial objects.

### Functions

- `1.gem2rds.ipynb`
  - Collects `.tissue.gem.gz` files from all samples.
  - Standardizes sample naming.
  - Batch-executes the conversion script.
- `stomics_seurat_only_gem2rds.R`
  - Reads the `.tissue.gem.gz` expression matrix.
  - Aggregates coordinates according to the specified `bin size`.
  - Builds a sparse expression matrix.
  - Creates a Seurat object with spatial coordinate information.
  - Outputs intermediate objects such as `*_bin30_seurat.rds`.

### Main Input

- Raw Stereo-seq files: `*.tissue.gem.gz`

### Main Output

- Seurat objects: `*_bin30_seurat.rds`

### Example Command

```bash
Rscript stomics_seurat_only_gem2rds.R \
  -i sample.tissue.gem.gz \
  -b 30 \
  -s GW6 \
  -o ./bin30Data
```

## Code 2: `2.selectCells.R` + `xselectCells.R`

This step is used to select the spatial regions or spots that will be retained for downstream analysis.

### Functions

- Reads the Seurat `rds` objects generated in Code 1.
- Stores `x` and `y` coordinates in a Seurat dimensional reduction object named `cod` for direct spatial visualization.
- Uses `FeaturePlot()` to inspect `nFeature_RNA` and marker gene expression across the tissue.
- Launches `xSelectCells()` for interactive manual region selection.
- Reads exported barcode files such as `Barcode_selected_*.csv`.
- Extracts the selected subset and saves it as a new Seurat `rds` object.

### Main Input

- `*_bin30_seurat.rds`

### Intermediate Output

- Manually selected barcode files: `Barcode_selected_*.csv`

### Main Output

- Filtered Seurat objects: `RDS/SeuratObj_*.rds`

### Notes

- `2.selectCells.R` contains example coordinate ranges and marker-based inspection logic for multiple developmental stages.
- `xselectCells.R` is an interactive selection utility that exports the selected barcodes.
- The current scripts still contain original cluster-specific paths and should be rewritten with relative paths or a configurable project root before public release.

## Code 2 (Optional): `2.(option)stereo_seq_spinal_cord_gw7.ipynb`

This is an optional step. It is intended for cases where marker genes and tissue morphology alone are insufficient to confidently define the target region. In such cases, NicheCompass can be used to perform unsupervised spatial clustering to identify biologically meaningful regions before manual or direct selection.

### Functions

- Loads the spatial transcriptomics object and prepares inputs for NicheCompass.
- Performs unsupervised clustering or niche partitioning on the spatial sample.
- Identifies regions that may correspond to biologically meaningful anatomical structures.
- Provides guidance for subsequent manual selection, or can be used for direct region selection based on clustering results.

### Recommended Use Cases

- Tissue boundaries are not visually clear.
- Marker gene signals are not sufficient to robustly define the target region.
- A data-driven spatial partition is preferred before deciding how to subset the sample.

### Main Input

- The spatial expression object generated in Code 1.

### Main Output

- Unsupervised spatial clustering or niche assignment results.
- Region annotations that can be used to guide or directly perform region selection.

### Notes

- The current directory includes a representative `GW7` notebook as an example of this optional step.
- This step is not required for all samples and should be used only when the target region cannot be confidently identified from prior knowledge alone.
- In practice, this step can be run first and then followed by `2.selectCells.R` for manual refinement, or the clustering results can be used directly to define the region of interest.

## Code 3: `3.convent.R`

This step converts the filtered Seurat objects into `h5ad` format for downstream Python-based analysis.

### Functions

- Scans all filtered Seurat objects under the `RDS/` directory.
- Converts each object to AnnData format using `sceasy::convertFormat()`.
- Writes the outputs to the `H5AD/` directory.

### Main Input

- `RDS/*.rds`

### Main Output

- `H5AD/*.h5ad`

## Recommended Execution Order

1. Prepare raw `.tissue.gem.gz` files for each sample.
2. Run `1.gem2rds.ipynb` and `stomics_seurat_only_gem2rds.R` to generate Seurat `rds` objects.
3. If the target region cannot be confidently defined using marker genes and morphology alone, run `2.(option)stereo_seq_spinal_cord_gw7.ipynb` to perform NicheCompass-based unsupervised spatial clustering.
4. Run `2.selectCells.R` together with `xselectCells.R` to manually select the target region, or directly select the region using the clustering results.
5. Verify that `Barcode_selected_*.csv` and the filtered `RDS/SeuratObj_*.rds` files are correct.
6. Run `3.convent.R` to batch-convert the filtered objects to `h5ad`.

## Outputs of This Section

After completing this section, two key types of intermediate files will be generated:

- Filtered Seurat objects for downstream R/Seurat-based analyses.
- Corresponding `h5ad` files for downstream Python/Scanpy or cross-platform workflows.
