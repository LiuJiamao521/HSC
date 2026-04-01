
export PATH=/cluster2/huanglab/jiamao/Apps/cellranger-8.0.1:$PATH

#GW7
cellranger count --id=GW7_gex \
                        --fastqs=/cluster/huanglab/lab/LabDataset/jiamao/HumanSpinalCord/RawData/scMulti-omics/GW7/RNA \
                        --sample=GW7 \
                        --transcriptome=/cluster2/huanglab/jiamao/Apps/refdata-gex-GRCh38-2024-A \
                        --create-bam=true \
                        --chemistry=ARC-v1 \
                        --localcores=8 \
                        --localmem=64