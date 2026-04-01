import os
import pandas as pd
import numpy as np
import scipy.sparse
from scipy import stats
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from upsetplot import from_contents, plot

# =========================
# 0. Global PDF export (editable vector text)
# =========================
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

def calculate_global_corr_worker(sub_df, gene_map, region_map, x_rna, x_atac):
    results = []
    for _, row in sub_df.iterrows():
        gene, region = row['target'], row['region']
        if gene not in gene_map or region not in region_map:
            continue
        g_idx, r_idx = gene_map[gene], region_map[region]
        g_expr = x_rna[:, g_idx]
        if scipy.sparse.issparse(g_expr): g_expr = g_expr.toarray().flatten()
        else: g_expr = np.ravel(g_expr)
        r_acc = x_atac[:, r_idx]
        if scipy.sparse.issparse(r_acc): r_acc = r_acc.toarray().flatten()
        else: r_acc = np.ravel(r_acc)
        
        if np.std(g_expr) == 0 or np.std(r_acc) == 0:
            rho, pval = 0.0, 1.0
        else:
            rho, pval = stats.spearmanr(r_acc, g_expr)
        results.append({'target': gene, 'region': region, 'global_rho': rho, 'global_pval': pval})
    return results

def run_lineage_comparison(
    adata_rna, 
    adata_atac, 
    region_df, 
    lineage_A, 
    lineage_B, 
    n_jobs=8,
    rho_threshold=0.03,
    pval_threshold=0.05,
    logfc_threshold=1.0
):
    # --- 1. Output paths ---
    name_a = lineage_A.replace('is_', '').replace('_lineage', '')
    name_b = lineage_B.replace('is_', '').replace('_lineage', '')
    comparison_name = f"{name_a}_vs_{name_b}"
    out_dir = os.path.join(".", "plot", comparison_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n>>> Running: {comparison_name}")

    # --- 2. Subset cells and compute global RNA–ATAC correlations ---
    mask = (adata_rna.obs[lineage_A] == 1) | (adata_rna.obs[lineage_B] == 1)
    sub_rna = adata_rna[mask].copy()
    sub_atac = adata_atac[mask].copy()
    X_rna = sub_rna.X.tocsc() if scipy.sparse.issparse(sub_rna.X) else sub_rna.X
    X_atac = sub_atac.X.tocsc() if scipy.sparse.issparse(sub_atac.X) else sub_atac.X
    gene_map = {name: i for i, name in enumerate(sub_rna.var_names)}
    region_map = {name: i for i, name in enumerate(sub_atac.var_names)}

    chunks = np.array_split(region_df, n_jobs * 4)
    results_list = Parallel(n_jobs=n_jobs)(
        delayed(calculate_global_corr_worker)(chunk, gene_map, region_map, X_rna, X_atac) for chunk in chunks
    )
    global_stats = pd.DataFrame([item for sublist in results_list for item in sublist])
    df_clean = pd.merge(region_df, global_stats, on=['target', 'region'], how='inner')
    links_validated = df_clean[(df_clean['global_rho'].abs() > rho_threshold) & (df_clean['global_pval'] < pval_threshold)].copy()

    # --- 3. Lineage-specific accessibility (log2 fold change) ---
    idx_A = np.where(sub_rna.obs[lineage_A] == 1)[0]
    idx_B = np.where(sub_rna.obs[lineage_B] == 1)[0]
    unique_regions = links_validated['region'].unique()
    region_indices = [region_map[r] for r in unique_regions if r in region_map]
    X_atac_sub = sub_atac.X[:, region_indices]
    mean_A = np.array(X_atac_sub[idx_A, :].mean(axis=0)).flatten()
    mean_B = np.array(X_atac_sub[idx_B, :].mean(axis=0)).flatten()
    log2fc = np.log2((mean_A + 1e-6) / (mean_B + 1e-6))
    
    df_abundance = pd.DataFrame({'region': [r for r in unique_regions if r in region_map], 'log2fc': log2fc})
    def classify_link(fc):
        if abs(fc) > logfc_threshold: return 'Specific'
        return 'Conserved'
    df_abundance['abundance_class'] = df_abundance['log2fc'].apply(classify_link)
    links_validated = pd.merge(links_validated, df_abundance, on='region', how='left')

    # --- 4. Donut / pie chart: Conserved vs Specific links ---
    counts = links_validated['abundance_class'].value_counts()
    n_cons, n_spec = counts.get('Conserved', 0), counts.get('Specific', 0)
    fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
    ax.pie([n_cons, n_spec], labels=None, autopct='%1.1f%%', startangle=90, colors=['#4A79A7', '#7DA87B'], 
           pctdistance=0.82, wedgeprops={'width': 0.5, 'edgecolor': 'white', 'linewidth': 2},
           textprops={'fontsize': 14, 'fontweight': 'bold', 'color': 'white'})
    ax.text(0, 0, f'{n_cons+n_spec:,}\nLinks', ha='center', va='center', fontsize=14, fontweight='bold')
    plt.title(f'Accessibility: {name_a} vs {name_b}', fontsize=15, fontweight='bold')
    plt.savefig(f"{out_dir}/{comparison_name}_Link_Pie.pdf", bbox_inches='tight', transparent=True)
    plt.close()

    # --- 5. Bar chart: per-gene category (replaces UpSet-style view) ---
    df_gene = links_validated.copy()
    df_gene["is_conserved"] = df_gene["abundance_class"] == "Conserved"
    df_gene["is_specific"] = df_gene["abundance_class"] == "Specific"

    gene_stats = df_gene.groupby("target").agg(
        total_links=('target', 'size'),
        cons_count=('is_conserved', 'sum'),
        spec_count=('is_specific', 'sum')
    ).reset_index()

    filtered_genes = gene_stats[gene_stats['total_links'] >= 10].copy()

    if not filtered_genes.empty:
        filtered_genes['pct_cons'] = filtered_genes['cons_count'] / filtered_genes['total_links']
        filtered_genes['pct_spec'] = filtered_genes['spec_count'] / filtered_genes['total_links']

        def classify_gene(row):
            if row['pct_cons'] >= 0.75: return 'Conserved Gene'
            elif row['pct_spec'] >= 0.75: return 'Specific Gene'
            else: return 'Dual Identity'

        filtered_genes['category'] = filtered_genes.apply(classify_gene, axis=1)

        cat_order = ['Dual Identity', 'Conserved Gene', 'Specific Gene']
        category_counts = filtered_genes['category'].value_counts().reindex(cat_order).fillna(0)

        plt.figure(figsize=(8, 6), dpi=300)
        colors = ['#9182C4', '#4A79A7', '#7DA87B']  # purple, blue, green
        
        ax = sns.barplot(x=category_counts.index, y=category_counts.values, palette=colors,
                         edgecolor='black', linewidth=1.2)
        
        for i, v in enumerate(category_counts.values):
            ax.text(i, v + (max(category_counts.values) * 0.02), str(int(v)), 
                    ha='center', fontweight='bold', fontsize=11)

        plt.title(f'Gene Classification: {name_a} vs {name_b}\n(Links >= 10, Threshold = 75%)', 
                  fontsize=14, fontweight='bold', pad=20)
        plt.ylabel('Number of Genes', fontsize=12, fontweight='bold')
        plt.xlabel('Regulatory Category', fontsize=12, fontweight='bold')
        sns.despine()
        plt.grid(axis='y', linestyle='--', alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(f"{out_dir}/{comparison_name}_Gene_Category_Bar.pdf", format='pdf', bbox_inches='tight')
        plt.close()
    else:
        print(f"Warning: {comparison_name}: no genes with >= 10 links; skipping gene-category bar plot.")

    links_validated.to_csv(f"{out_dir}/{comparison_name}_results.csv", index=False)
    print(f">>> Finished: {comparison_name}")
