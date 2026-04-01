# Preprocessing: Cell Ranger ARC and library configuration

Command-line examples for **10x Multiome (GEX + ATAC)** quantification from raw FASTQs. No downstream Python environment is required for this step.

## Workflow

1. **`1.cellranger-arc.sh`**  
   Runs **Cell Ranger ARC** `count` on FASTQs to produce joint **GEX + ATAC** `outs/` (e.g. `filtered_feature_bc_matrix`, `atac_fragments.tsv.gz`). Edit `PATH`, `--reference`, `--libraries`, `--id`, `--localcores`, and `--localmem` for your cluster.

2. **`GW7.csv`**  
   Example **libraries CSV** for Cell Ranger ARC (GEX/ATAC library paths and sample names), referenced by `--libraries`.

## Input and output

- **Input:** Raw FASTQs, ARC reference package, libraries table CSV.  
- **Output:** Per-sample `outs/` from `cellranger-arc count` (matrices, fragments, QC HTML, etc.).

