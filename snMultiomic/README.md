# Single-nucleus multiome (snMultiomic) workflows

This directory contains the **single-nucleus RNA + ATAC (10x Multiome)** analysis framework for this project: raw quantification with Cell Ranger ARC, Scanpy and SnapATAC2 for modality-specific processing, motif enrichment, RNA velocity, multiscale ATAC footprinting, SCENIC+/pyCistopic-based gene regulatory networks, and CellOracle in silico TF perturbation.

Each numbered submodule includes its own **`README.md`** with script- and notebook-level details (same pattern as `spatial/README.md`).

## Module overview

1. **`01_preprocessing_cellranger/`**  
   Shell templates for **Cell Ranger ARC** joint GEX + ATAC counting and an example **`GW7.csv`** libraries file; primary outputs are `filtered_feature_bc_matrix` and `atac_fragments.tsv.gz` under each sample `outs/`.

2. **`02_snrna_scanpy/`**  
   **Scanpy** notebooks: load and merge multi-sample GEX matrices, QC and embedding (Part 1), integrate metadata and cell-type labels (Part 2), differential expression and marker visualization (Part 3).

3. **`03_snatac_snapatac2/`**  
   **SnapATAC2** notebooks for multi-sample fragment import, peak matrices and integration, **label transfer** from snRNA, and peak accessibility summaries.

4. **`04_motif_homer/`**  
   **HOMER**-oriented notebooks: build CRE **BED** sets by cell type, **shuffle** backgrounds, run/compare motif enrichment with **module-to-subgroup** summaries.

5. **`05_rna_velocity_scvelo/`**  
   **Cell Ranger** (BAM-enabled) and **velocyto** shell steps, then **scVelo** notebooks for dynamics, latent time, and an optional oligo practice notebook.

6. **`06_footprint_scprinter_seq2print/`**  
   **scPrinter** workflows for Tn5 insertions, TF binding scores, **multiscale footprints**, plus **Seq2Print** attribution visualization against peaks.

7. **`07_grn_scenicplus_pycistopic/`**  
   **pyCistopic** topic modeling from peak matrices, **BigWig** generation by lineage, full **SCENIC+** pipeline (`h5mu` outputs), gene–peak link plots, **eGRN** network export, and **`lineage_analysis.py`** for lineage-wise CRE–gene validation.

8. **`08_celloracle_perturbation/`**  
   **CellOracle** notebooks: peaks to TF motif / **Base-GRN**, **GRN** fitting and network analysis from scRNA + Base-GRN, **in silico TF perturbation** and vector-field visualization.

## Recommended end-to-end order

1. Quantify raw Multiome data (`01`).
2. Run snRNA analysis (`02`) and snATAC analysis (`03`) in parallel where practical; run **label transfer** in `03` once RNA annotations are stable.
3. Prepare peaks and run motif enrichment (`04`).
4. If using velocity, generate BAMs and loom, then run scVelo (`05`)—can overlap with `02`/`03` once matrices exist.
5. Optional deep ATAC interpretation: footprinting (`06`) and GRN inference (`07`)—both benefit from clean fragments/peaks and consistent RNA metadata.
6. Run CellOracle perturbation (`08`) after RNA preprocessing and Base-GRN–related inputs are available.

## Data products across the snMultiomic pipeline

- **Core objects:** AnnData **`h5ad`** (RNA and ATAC), SnapATAC2 peak matrices, SCENIC+ **`h5mu`** / MuData-style multi-modal objects.
- **Intermediate files:** Cell Ranger `outs/`, velocyto **loom**, per-lineage **BigWig** tracks, HOMER foreground/background BED and results tables.
- **Regulatory and dynamics outputs:** HOMER motif summaries, **eRegulons** and gene–peak links, multiscale footprint and binding-score artifacts, CellOracle **Links** / simulation vector plots, scVelo latent-time and stream figures.

## Notes

- Paths, sample IDs, and tool versions embedded in notebooks are **environment-specific**; replace them before reuse.
- For software stacks, favor dedicated conda/env per heavy module (e.g. SCENIC+, CellOracle, scPrinter) as implied by the original notebooks.
