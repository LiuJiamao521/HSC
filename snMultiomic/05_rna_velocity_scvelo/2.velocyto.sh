cd /cluster/huanglab/lab/LabDataset/jiamao/HumanSpinalCord/ProcessData/Cellranger/

velocyto run -b /cluster/huanglab/lab/LabDataset/dhong/HumanSpinalCord/ProcessData/Cellranger/GW7/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o /cluster/huanglab/lab/LabDataset/jiamao/HumanSpinalCord/ProcessData/Cellranger/GW7/ -m hg38_rmsk.gtf /cluster/huanglab/lab/LabDataset/dhong/HumanSpinalCord/ProcessData/Cellranger/GW7/outs/gex_possorted_bam.bam genes.gtf
