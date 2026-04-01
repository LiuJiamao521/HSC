# In silico perturbation: CellOracle Base-GRN, networks, and TF KO simulation

Follows **CellOracle** tutorials: scan **TF motifs** on **scATAC peaks** to build a **Base-GRN**, fit **cluster-specific GRNs** with **scRNA**, then run **in silico TF knockdown/overexpression** and visualize transition vectors.

## Workflow

1. **`01_atac_peaks_to_TFinfo.ipynb`**  
   Motif scanning on ATAC peaks to produce CellOracle **TF info / Base-GRN** inputs (see CellOracle *atac_peaks_to_TFinfo* documentation).

2. **`02_network_analysis.ipynb`**  
   Fit **per-cluster GRNs** (**Links** objects) from RNA + Base-GRN; compare network structure and export for external graph tools.

3. **`03_tf_ko_simulation.ipynb`**  
   Load Oracle objects and inferred GRNs, fit models for simulation, run **TF perturbations**, compute **transition vectors**, plot **quiver / vector fields**, optionally compare to developmental directions (see notebook sections).

## Input and output

- **Input:** Preprocessed **expression** in CellOracle-compatible form, peaks, motif databases, species-appropriate Base-GRN settings.  
- **Output:** Base-GRN, Links/Oracle objects, perturbation vector plots and intermediate tables.
