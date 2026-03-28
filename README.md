# Human Spinal Cord Multi-omics Atlas Code Release

## Overview
This repository organizes the analysis code used to study the developing human spinal cord through snMultiome and Stereo-seq profiling. The release is structured to match the manuscript workflow and to make it easier for reviewers and readers to trace each major analysis module, from preprocessing and lineage analysis to spatial mapping, regulatory network inference, niche definition, and trait enrichment.

## Repository Structure
```text
Project_Code_Release/
├── snMultiomic/
│   ├── 01_preprocessing_cellranger_velocyto/
│   ├── 02_scrna_de_go/
│   ├── 03_atac_peaks_cCRE/
│   ├── 04_dars_snapatac2/
│   ├── 05_motif_homer/
│   ├── 06_rna_velocity_scvelo/
│   ├── 07_footprinting_scprinter_seq2print/
│   ├── 08_grn_scenicplus_pycistopic/
│   └── 09_celloracle_perturbation/
├── spatial/
│   ├── 01_stereo_preprocessing/
│   ├── 02_integration_seurat_harmony_knn/
│   ├── 03_cell2location/
│   ├── 04_spatial_multimapping_atac_magic_chromvar/
│   ├── 05_nichecompass_niches/
│   ├── 06_neighborhood_ncem_cellchat/
│   └── 07_trait_gwas_gsmap_ldsc/
├── plotting/
│   ├── main_figures/
│   └── extended_data/
├── shared/
│   └── resources/
├── method.md
├── prompt.txt
├── .gitignore
├── LICENSE
└── README.md
```

## Workflow
1. Process snMultiome RNA/ATAC data and perform cell annotation.
2. Identify marker genes, cCREs, DARs, enriched motifs, and lineage dynamics.
3. Infer enhancer-driven gene regulatory networks and simulate TF perturbation.
4. Process Stereo-seq data and integrate it with stage-matched single-cell references.
5. Map cell types to tissue space with Seurat/Harmony/kNN and cell2location.
6. Reconstruct spatial TF activity, regulatory programs, niches, and cell-cell communication.
7. Quantify spatial and cell-type-level trait enrichment and generate figure panels.

## Reproducibility
Representative project notebooks and scripts were copied from the original `ATAC/` and `ST/` analysis trees into logic-based folders. The release keeps code and workflow organization, but intentionally excludes large raw or intermediate datasets such as `.h5ad`, `.rds`, `.fastq`, fragment files, and big matrices.

To reproduce the analyses:
1. Prepare the required input objects in a local `data/` directory outside this repository.
2. Review the relevant module under `snMultiomic/` or `spatial/`.
3. Update input and output paths in the copied notebooks/scripts to your environment.
4. Run the module notebooks/scripts in workflow order.
5. Use the scripts in `plotting/` to regenerate manuscript figures from processed intermediate results.

## Path Conventions
The original analysis notebooks were developed on an internal cluster and may still contain absolute paths such as `/cluster2/huanglab/...`. These paths should be replaced with your own project root, environment variables, or relative paths before public release or reuse. A recommended pattern is to define a project root once at the top of each script or notebook and build all file paths from that root.

## Notes
- `snMultiomic/01_preprocessing_cellranger_velocyto/` contains downstream notebooks representing the preprocessing stage; dedicated Cell Ranger ARC and velocyto submission scripts were not found in the currently available source tree.
- `spatial/07_trait_gwas_gsmap_ldsc/` currently includes the available LDSC-related code from the source tree. The manuscript also describes gsMap-based spatial trait mapping; if those scripts are stored elsewhere, they should be added in a later update.
- `shared/resources/copied_files_inventory.txt` records the files copied into this release.

## Citation
Citation information will be added upon manuscript publication.
