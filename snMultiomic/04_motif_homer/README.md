# Motif: candidate CREs and HOMER enrichment

Prepare **BED** files for **cell-type- or contrast-specific candidate regulatory elements (CREs / peaks)** and run **HOMER** motif enrichment with **module-to-subgroup** summaries; includes **shuffle** backgrounds to assess specificity.

## Workflow

1. **`1.cres_bed.ipynb`**  
   Aggregate **CRE ↔ cell_type** from peaks and labels; export **BED** or grouped region files for HOMER.

2. **`2.shuffle_background.ipynb`**  
   Build **shuffle / control** region sets for comparison with foreground CREs to reduce composition-driven false positives.

3. **`3.homer_module2subgroup.ipynb`**  
   Organize and visualize HOMER results from **modules** to **subgroups** to interpret which motif programs enrich in which populations.

## Input and output

- **Input:** Annotated peak sets, cell grouping, reference genome (for HOMER).  
- **Output:** BED files, HOMER statistics, summary tables and plots.
