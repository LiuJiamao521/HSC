#!/usr/bin/env bash
# Cell Ranger ARC: joint Gene Expression + ATAC from FASTQs.
# Edit PATH, --reference, --libraries, --id, and resource flags for your environment.

export PATH=/cluster2/huanglab/jiamao/Apps/cellranger-arc-2.0.2:$PATH

cellranger-arc count --id=GW7 \
                       --reference=/cluster2/huanglab/jiamao/Apps/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/cluster/huanglab/lab/LabDataset/jiamao/HumanSpinalCord/ProcessData/Cellranger/GW7.csv \
                       --localcores=8 \
                       --localmem=64
