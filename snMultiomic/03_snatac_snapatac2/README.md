# snATAC: SnapATAC2 multi-sample pipeline

**SnapATAC2** processing of **Cell Ranger ARC** `fragments`: import, QC, peak calling / count matrix, embedding and clustering, **cross-modal label transfer** from snRNA, and **peak-level accessibility** analysis.

## Workflow

1. **`1.snapATAC2.ipynb`**  
   Multi-sample human spinal cord snATAC pipeline: load `atac_fragments.tsv.gz`, QC, peaks / matrix, integration and visualization (*Multi-sample Pipeline: analyzing snATAC-seq data of human spine samples*).

2. **`2.snATAC_annotation.ipynb`**  
   Annotate ATAC cells using snRNA-derived types; filtering, UMAP, marker visualization in ATAC space, sample distribution.

3. **`3.Peak_accessbility.ipynb`**  
   **Peak / chromatin accessibility** summaries and plots on top of annotations and peak matrices (often by cell type or contrast).

## Input and output

- **Input:** Per-sample `outs/atac_fragments.tsv.gz` (optional annotated RNA `h5ad` for label transfer).  
- **Output:** SnapATAC2 `AnnData` (peak × cell), intermediates, figures; feeds motif, footprint, pyCistopic / SCENIC+.
