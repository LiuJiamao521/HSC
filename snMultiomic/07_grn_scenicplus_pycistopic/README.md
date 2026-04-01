# Gene regulatory networks: pyCistopic, SCENIC+, and eGRN

From **peak × cell** matrices to **cis-regulatory topics**, **enhancer–gene links**, and the full **SCENIC+** pipeline, plus **eGRN** network construction/visualization and **lineage-wise CRE–gene link** statistics (`lineage_analysis.py`).

## Workflow

1. **`1.pycistopic.ipynb`**  
   Build **cisTopic** from SnapATAC2 **peak matrix**, run LDA (e.g. MALLET), obtain topics and cell–topic distributions for downstream GRN steps.

2. **`2.bigwig_generate.ipynb`**  
   Aggregate **BigWigs** (or track files) by **lineage** or group from fragments for SCENIC+, genome browsers, and integration.

3. **`3.scenicplus_workflow.ipynb`**  
   Configure and run **SCENIC+**; main artifact **`scplusmdata.h5mu`** (eRegulons and enrichment scores). Includes eRegulon embedding, specificity scores, dot heatmaps, and **track plots**.

4. **`4.cre_gene_links.ipynb`** / **`5.links_compare.ipynb`**  
   Plot and compare **gene–peak links** along trajectories or by batch/group to sanity-check regulatory inference.

5. **`6.egrn_construct.ipynb`**  
   Build and visualize **eGRNs** from SCENIC+ outputs with **networkx** (export to Cytoscape).

6. **`lineage_analysis.py`**  
   For two **lineages**, compute global **Spearman** correlation between CRE–gene links on matched RNA/ATAC matrices, combine with **lineage-specific accessibility log2FC**, filter, and plot (parallel workers, Upset-style outputs).

## Input and output

- **Input:** peak `h5ad`, fragments, RNA `h5ad`, SCENIC+ config and resources (motifs, rankings, etc.).  
- **Output:** cisTopic runs, BigWigs, `h5mu`, eRegulon tables, network figures, lineage comparison stats and PDFs.

## Notes

- pyCistopic/MALLET and SCENIC+ are resource-heavy; localize paths and versions in notebooks.  
- `lineage_analysis.py` expects boolean lineage columns in `adata.obs` (e.g. `is_*_lineage`).
