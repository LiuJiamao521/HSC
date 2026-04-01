#!/bin/bash
# Cell-type specific CREs
python run_LDSC_nat.py \
  -i /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/bed \
  -o /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/output_nat \
  -q /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/qsub_nat \
  -c /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/config.txt \
  -r

# Script for summary 
input_dir="/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/output_nat/04_Enrichment_result"
output_file="hsc.summary.tsv"
> "$output_file"
for file in "$input_dir"/*.cell_type_results.txt; do
    sample=$(basename "$file" | cut -d'.' -f1)
    awk -v sample="$sample" 'NR!=1 {print $0 "\t" sample}' "$file"
done >> "$output_file"


# # Cell-type specific CREs
# for bed_file in /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/bed/*.bed; do
#   filename=$(basename "$bed_file" .bed)
#   python run_LDSC_cell.py \
#     -i /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/bed \
#     -p "$filename" \
#     -o /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/output2 \
#     -q /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/qsub2 \
#     -c /cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/GWAS/config.txt
# done