# Footprinting: scPrinter multiscale footprints and Seq2Print attributions

**scPrinter (`scp`)** for **Tn5 insertion** profiles at single-cell or pseudobulk resolution, **TF binding site scoring**, and **multiscale footprint** metrics; some notebooks integrate **Seq2Print** **attribution** `.npz` outputs for sequence-level interpretation.

## Workflow

1. **`1.HSC_snATAC_oligo.ipynb`**  
   Example **HSC / oligolineage** snATAC: import fragments, group barcodes, basic scPrinter workflow (insertion tracks, group comparisons).

2. **`2.footprint_oligo.ipynb`**  
   **Multiscale footprint with scprinter:** fragments → insertions → load dispersion / binding-score models → **TFBS scores** → **multiscale footprints** and synced visualization (includes bigWig prep notes).

3. **`3.motif_oligo.ipynb`**  
   Load Seq2Print-style **`.npz` attributions** and **peak BED**; **sequence attribution** plots in motif/region context (oligo-focused analysis line).

## Input and output

- **Input:** `fragments.tsv.gz`, barcode groups, pretrained binding/footprint model paths; optional Seq2Print attributions and region BED.  
- **Output:** Insertion and footprint plots, binding-score and multiscale matrix-like objects, publication figures.

