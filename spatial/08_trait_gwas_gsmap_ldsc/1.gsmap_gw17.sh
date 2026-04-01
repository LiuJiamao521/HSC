#!/bin/bash

# data preprocessing, just run once
gsmap run_find_latent_representations \
    --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
    --sample_name 'GW17' \
    --input_hdf5_path '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/gsMap_example_data/ST/GW17.h5ad' \
    --data_layer 'counts' \
    --annotation 'annotation'

gsmap run_latent_to_gene \
    --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
    --sample_name 'GW17' \
    --num_neighbour 30 \
    --num_neighbour_spatial 50 \
    --latent_representation 'latent_GVAE' \
    --annotation 'annotation'

for CHROM in {1..22}
do
    gsmap run_generate_ldscore \
        --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
        --sample_name 'GW17' \
        --chrom $CHROM \
        --bfile_root 'gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink_hg38/1000G.EUR.QC' \
        --gtf_annotation_file 'gsMap_resource/genome_annotation/gtf/gencode.v46.basic.annotation.gtf' \
        --gene_window_size 100000 \
        --enhancer_annotation_file 'gsMap_resource/genome_annotation/enhancer/by_tissue_hg38/ALL/ABC_roadmap_merged.bed' \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'enhancer_first' \
        --spots_per_chunk 300
done

# gsmap run_spatial_ldsc \
#     --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
#     --sample_name 'GW17' \
#     --trait_name 'SCZ' \
#     --sumstats_file '/cluster/huanglab/lab/LabDataset/jiamao/GWAS/SCZ/SCZ.sumstats.gz' \
#     --num_processes 8
# #    --w_file 'gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights.' \


# gsmap run_cauchy_combination \
#     --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
#     --sample_name 'GW17' \
#     --trait_name 'SCZ' \
#     --annotation 'annotation'

# gsmap run_report \
#     --workdir '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine' \
#     --sample_name 'GW17' \
#     --trait_name 'SCZ' \
#     --annotation 'annotation' \
#     --sumstats_file '/cluster/huanglab/lab/LabDataset/jiamao/GWAS/SCZ/SCZ.sumstats.gz' \
#     --top_corr_genes 50

# run the script for multiple traits using for loop

# Set working directory and GWAS data directory
WORKDIR="/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/gsMap/Human_Spine"
GWAS_DIR="/cluster3/labData/jiamao/GWAS"
SAMPLE_NAME="GW17"
NUM_PROCESSES=8
ANNOTATION="annotation"

# Define trait list to process
TRAIT_LIST=("AHD" "AR" "ATM" "Benign" "CRPS" "Encephalitis" "Epilepsy" "GGE" "Insomnia" "Parasomnia" "Depression" "BP" "Neuroticism" "Height" "IQ" "SPA" "AD" "AN" "ASD" "MS" "PTSD" "ADHD" "ANX" "SCZ" "CDG" "MDD" "OCD" "ALS")

# Iterate through trait list
for trait_name in "${TRAIT_LIST[@]}"; do
    # Build sumstats file path
    sumstats_file="$GWAS_DIR/$trait_name/${trait_name}.sumstats.gz"
    
    # Check whether the sumstats file exists
    if [ ! -f "$sumstats_file" ]; then
        echo "Warning: Sumstats file not found for trait $trait_name at $sumstats_file, skipping..."
        continue
    fi
    
    echo "Processing trait: $trait_name"
    
    # Run run_spatial_ldsc command
    if ! gsmap run_spatial_ldsc \
        --workdir "$WORKDIR" \
        --sample_name "$SAMPLE_NAME" \
        --trait_name "$trait_name" \
        --sumstats_file "$sumstats_file" \
        --num_processes "$NUM_PROCESSES"; then
        echo "Error occurred in run_spatial_ldsc for trait $trait_name, skipping..."
        continue
    fi
    
    # Run run_cauchy_combination command
    if ! gsmap run_cauchy_combination \
        --workdir "$WORKDIR" \
        --sample_name "$SAMPLE_NAME" \
        --trait_name "$trait_name" \
        --annotation "$ANNOTATION"; then
        echo "Error occurred in run_cauchy_combination for trait $trait_name, skipping..."
        continue
    fi
    
    # Run run_report command
    if ! gsmap run_report \
        --workdir "$WORKDIR" \
        --sample_name "$SAMPLE_NAME" \
        --trait_name "$trait_name" \
        --annotation "$ANNOTATION" \
        --sumstats_file "$sumstats_file" \
        --top_corr_genes 50; then
        echo "Error occurred in run_report for trait $trait_name, skipping..."
        continue
    fi
    
    echo "Successfully processed trait: $trait_name"
done

echo "All specified traits processed."
