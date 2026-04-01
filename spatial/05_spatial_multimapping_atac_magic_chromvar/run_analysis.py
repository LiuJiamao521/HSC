import os
os.environ["SCPRINTER_DATA"] = "/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ATAC/scPrinter/scPrinter/scprinter/docs/reference/Datasets"
import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import pickle
import torch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from scanpy.plotting.palettes import zeileis_28
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map, thread_map
import anndata
import scanpy as sc
import statistics as stat
import json
import csv
import re
import copy
from sklearn.preprocessing import OneHotEncoder

# We can calculate chromVAR motif scores using either GPU (device = "cuda", much faster) or CPU (device = "cpu", slower)
device = "cuda:0"

if device == "cuda:0":
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc
    import anndata
    import cupy as cp
    import cupyx as cpx
    import time
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True, # Allows oversubscription
        pool_allocator=True, # default is False
        devices=0, # GPU device IDs to register. By default registers only GPU 0.
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

# Set paths and base configuration
data_dir = '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/visualization3/reconstructed_SpATAC'
work_dir = '/cluster2/huanglab/jiamao/Project/HumanSpinalCord/Work/P3/ST/chromvar'
genome = scp.genome.hg38

# Create output directory if it does not exist
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

# 2. Collect all .h5ad files in the directory
files = [f for f in os.listdir(data_dir) if f.endswith('.h5ad') and f.startswith('Spatial_ATAC_')][3:4]
print(f"Found files: {files}")

# 3. Start iterative processing
for file_name in files:
    # Extract sample ID (e.g., 'Spatial_ATAC_GW6.h5ad' -> 'GW6')
    sample_id = file_name.replace('Spatial_ATAC_', '').replace('.h5ad', '')
    print(f"\n--- Processing Sample: {sample_id} ---")
    
    try:
        # Read data
        adata_path = os.path.join(data_dir, file_name)
        adata = sc.read_h5ad(adata_path)
        
        # Filter peaks with zero coverage
        coverage = adata.X.sum(axis=0)
        adata = adata[:, coverage > 0].copy()

        # Sample background peaks
        scp.chromvar.sample_bg_peaks(adata,
                                     genome=genome,
                                     method='chromvar',
                                     niterations=250)

        # Scan motifs
        # Note: tune n_jobs based on available CPU cores
        motif = scp.motifs.FigR_Human_Motifs(genome,
                                             bg=list(adata.uns['bg_freq']),
                                             n_jobs=100,
                                             pvalue=5e-5, 
                                             mode='motifmatchr')
        motif.prep_scanner(None, pvalue=5e-5)
        motif.chromvar_scan(adata)

        # Compute motif scores
        # chunk_size can be tuned for GPU/CPU memory capacity
        chromvar = scp.chromvar.compute_deviations(adata, chunk_size=5000, device=device)
        
        # Transfer spatial coordinates
        if 'spatial' in adata.obsm:
            chromvar.obsm['spatial'] = adata.obsm['spatial']
        
        # Save result
        output_path = os.path.join(work_dir, f'chromvar_cisbp_{sample_id}.h5ad')
        chromvar.write(output_path)
        print(f"Successfully saved: {output_path}")

    except Exception as e:
        print(f"Error processing {sample_id}: {str(e)}")
        continue

print("\n--- All samples finished ---")