# RNA velocity: Cell Ranger BAM, velocyto, and scVelo

End-to-end **RNA velocity**: **spliced / unspliced** counting from BAMs and **scVelo** dynamics—Cell Ranger with BAM → loom → `AnnData` and latent-time / stream plots.

## Workflow

1. **`1.cellranger.sh`**  
   **Cell Ranger count** on GEX FASTQs with **`--create-bam=true`** so velocyto can count intronic/exonic reads (match `chemistry` etc. to your ARC library).

2. **`2.velocyto.sh`**  
   **velocyto run** on BAM + gene annotation + repeat mask to produce **loom** (or compatible output); align output paths with your pipeline.

3. **`3.scvelo_ldata.ipynb`**  
   Load annotated `h5ad` with velocity layers; **scVelo** for dynamics, latent time, stream / phase plots.

4. **`4.practice_oligo.ipynb`**  
   **Practice / subset** notebook for **oligo** (or similar) to tune parameters separately from the main flow.

## Input and output

- **Input:** GEX BAM, GTF, repeat mask (velocyto); annotated `h5ad` (scVelo).  
- **Output:** loom, objects with `layers['spliced'/'unspliced']`, velocity figures.

