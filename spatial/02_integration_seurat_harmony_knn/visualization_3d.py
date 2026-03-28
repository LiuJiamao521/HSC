# visualization_3d.py
"""
3D visualization utility for multi-omics integration.

Dependencies:
- numpy
- pandas
- matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def visualization_3d(
    rna_umap: pd.DataFrame,
    rna_anno: pd.DataFrame,
    spatial_coord: pd.DataFrame,
    link_st_matrix: pd.DataFrame,
    link_atac_matrix: pd.DataFrame,
    atac_umap: pd.DataFrame,
    color_map: dict,
    target_subtypes: list,
    # Coordinate processing parameters
    rna_rot: int = 0,
    spatial_rot: int = 0,
    atac_rot: int = 0,
    scale_rna: float = 1.0,
    scale_spatial: float = 0.6,
    scale_atac: float = 1.0,
    scale_spot: float = 1.0,
    # Visualization parameters
    figsize: tuple = (16, 12),
    elev: float = 30,
    azim: float = 300,
    title: str = "Spatial-GEX-ATAC Integration",
    line_alpha: float = 0.3,
    point_sizes: dict = {'rna': 5, 'atac': 5, 'spatial': 25},
    axis_scale: float = 0.08
):
    """
    Visualize the integration of spatial transcriptomics, RNA-seq, and ATAC-seq in 3D.

    Args:
        rna_umap (pd.DataFrame):
            RNA UMAP coordinates indexed by barcode, containing columns
            ['UMAP_1', 'UMAP_2'].
        rna_anno (pd.DataFrame):
            RNA cell-type annotations indexed consistently with `rna_umap`,
            containing a ['cell_name'] column.
        spatial_coord (pd.DataFrame):
            Spatial coordinates indexed by barcode, containing columns
            ['cod1', 'cod2'].
        link_st_matrix (pd.DataFrame):
            Spatial-RNA linkage table. The index is not important, but the table
            should contain:
            - id_st: spatial spot ID
            - id_rna: RNA cell ID
            - cluster_rna: cell type
        link_atac_matrix (pd.DataFrame):
            ATAC-RNA linkage table. The index is not important, but the table
            should contain:
            - id_atac: ATAC cell ID
            - id_rna: RNA cell ID
            - cluster_rna: cell type
        atac_umap (pd.DataFrame):
            ATAC UMAP coordinates indexed by barcode, containing columns
            ['UMAP_1', 'UMAP_2'].
        color_map (dict):
            Mapping from cell type to color.
        target_subtypes (list):
            List of highlighted cell types to display.

        rna_rot (int, optional):
            Rotation angle for RNA coordinates. Default is 0 degrees and usually
            no rotation is needed.
        spatial_rot (int, optional):
            Rotation angle for spatial coordinates. This should be adjusted based
            on the dorsoventral orientation of the tissue section.
        atac_rot (int, optional):
            Rotation angle for ATAC coordinates. Default is 0 degrees and usually
            no rotation is needed.
        scale_rna (float, optional):
            Scaling factor for RNA coordinates. Default is 1.0.
        scale_spatial (float, optional):
            Scaling factor for spatial coordinates. Default is 0.6, mainly to
            shrink the tissue section for clearer display.
        scale_atac (float, optional):
            Scaling factor for ATAC coordinates. Default is 1.0.
        scale_spot (float, optional):
            Scaling factor for spatial point size. Default is 1.0.

        figsize (tuple, optional):
            Figure size. Default is (16, 12).
        elev (float, optional):
            Elevation angle of the 3D view. Default is 30 degrees.
        azim (float, optional):
            Azimuth angle of the 3D view. Default is 300 degrees.
        title (str, optional):
            Figure title. Default is "Spatial-GEX-ATAC Integration".
        line_alpha (float, optional):
            Line transparency. Default is 0.3.
        point_sizes (dict, optional):
            Point sizes for each layer. Default is:
            {'rna': 5, 'atac': 5, 'spatial': 25}
        axis_scale (float, optional):
            Axis length scale. Default is 0.08.

    Returns:
        matplotlib.figure.Figure:
            Figure object.
    """
    
    # ================== Coordinate processing ==================
    def safe_normalize(arr, scale=1.0):
        centroid = np.mean(arr, axis=0)
        centered = arr - centroid
        max_abs = np.max(np.abs(centered), axis=0)
        max_abs = np.where(max_abs < 1e-8, 1.0, max_abs)
        return (centered * (0.5 / max_abs)) * scale

    def rotate_coords(coords, angle):
        theta = np.radians(angle)
        rot_matrix = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ])
        return coords @ rot_matrix

    # Process coordinates for each layer
    def process_layer(data, cods, rot, scale):
        base = data[cods].values
        centroid = np.mean(base, axis=0)
        rotated = rotate_coords(base - centroid, rot) + centroid
        return safe_normalize(rotated, scale)

    # Process RNA coordinates
    umap_xy = process_layer(rna_umap, ['UMAP_1', 'UMAP_2'], rna_rot, scale_rna)
    
    # Process spatial coordinates
    spatial_xy = process_layer(spatial_coord, ['cod1', 'cod2'], spatial_rot, scale_spatial)
    
    # Process ATAC coordinates
    atac_xy = process_layer(atac_umap, ['UMAP_1', 'UMAP_2'], atac_rot, scale_atac)

    # ================== Data preparation ==================
    spatial_clusters = link_st_matrix.set_index('id_st')['cluster_rna']
    atac_clusters = link_atac_matrix.set_index('id_atac')['cluster_rna']
    target_color_map = {key: color_map[key] for key in target_subtypes if key in color_map}

    # Generate the filtered set of linked spatial IDs
    valid_link_st = link_st_matrix[link_st_matrix['cluster_rna'].isin(target_subtypes)]['id_st'].unique()
    valid_link_st = set(valid_link_st).intersection(spatial_coord.index)  # Keep valid IDs only
    
    # Generate dynamic point-size array
    spatial_sizes = [
        point_sizes['spatial'] * scale_spot if idx in valid_link_st  # Enlarge points that satisfy the condition
        else point_sizes['spatial']
        for idx in spatial_coord.index
    ]

    # ================== Visualization settings ==================
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    ax.computed_zorder = False
    
    # Set coordinate ranges
    ax.set_zlim(-1.2, 1.2)
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)

    # ================== Draw layers ==================
    # ATAC layer (z=-1)
    ax.scatter(
        atac_xy[:,0], atac_xy[:,1], -1,
        c=[target_color_map.get(atac_clusters.get(x, None), '#D3D3D3') for x in atac_umap.index],
        s=point_sizes['atac'], alpha=0.7, edgecolor='w', linewidth=0.3,
        zorder=1
    )

    # ATAC-RNA links
    for _, row in link_atac_matrix.iterrows():
        try:
            atac_idx = atac_umap.index.get_loc(row['id_atac'])
            rna_idx = rna_umap.index.get_loc(row['id_rna'])
            color = target_color_map.get(row['cluster_rna'], '#D3D3D3')
            if color != '#D3D3D3':
                ax.plot(*zip(atac_xy[atac_idx], umap_xy[rna_idx]), [-1, 0],
                        color=color, linewidth=0.5, alpha=line_alpha, zorder=2)
        except KeyError:
            continue

    # RNA layer (z=0)
    ax.scatter(
        umap_xy[:,0], umap_xy[:,1], 0,
        c=[target_color_map.get(ct, '#D3D3D3') for ct in rna_anno],
        s=point_sizes['rna'], alpha=0.7, edgecolor='w', linewidth=0.3,
        zorder=3
    )

    # Spatial-RNA links
    for _, row in link_st_matrix.iterrows():
        try:
            st_idx = spatial_coord.index.get_loc(row['id_st'])
            rna_idx = rna_umap.index.get_loc(row['id_rna'])
            color = target_color_map.get(row['cluster_rna'], '#D3D3D3')
            if color != '#D3D3D3':
                ax.plot(*zip(spatial_xy[st_idx], umap_xy[rna_idx]), [1, 0],
                        color=color, linewidth=0.5, alpha=line_alpha, zorder=4)
        except KeyError:
            continue

    # Spatial layer (z=1)
    ax.scatter(
        spatial_xy[:,0], spatial_xy[:,1], 1,
        c=[target_color_map.get(spatial_clusters.get(x, None), '#D3D3D3') for x in spatial_coord.index],
        s=spatial_sizes, alpha=0.8, edgecolor='w', linewidth=0.5,
        zorder=5
    )

    # ================== Axis system ==================
    axis_style = {
        'color': '#404040', 
        'linewidth': 1.2, 
        'arrow_length_ratio': 0.15,
        'alpha': 0.9
    }

    label_style = {
        'fontsize': 9,
        'fontweight': 'bold',
        'bbox': dict(boxstyle="round,pad=0.2", facecolor='white', edgecolor='none', alpha=0.9)
    }

    # Draw layer-specific axes
    def draw_axes(origin, z, labels, axis_length=0.08):
        # UMAP axes
        ax.quiver(origin[0], origin[1], z, axis_length, 0, 0, **axis_style)
        ax.quiver(origin[0], origin[1], z, 0, axis_length, 0, **axis_style)
        # Labels
        ax.text(origin[0]+axis_length*1.2, origin[1], z, labels[0],
                color=axis_style['color'], rotation=-15, ha='left', va='center', **label_style)
        ax.text(origin[0], origin[1]+axis_length*1.2, z, labels[1], 
                color=axis_style['color'], rotation=25, ha='center', va='bottom', **label_style)

    # RNA-layer axes
    draw_axes([-0.48, -0.48], 0, ['UMAP1', 'UMAP2'], axis_scale)
    ax.text(0.5, 0, 0, 'snRNA-seq', **label_style)

    # Spatial-layer axes
    origin = [-0.35, -0.35]
    ax.quiver(origin[0], origin[1], 1, -axis_scale, 0, 0, **axis_style)
    ax.quiver(origin[0], origin[1], 1, axis_scale, 0, 0, **axis_style)
    ax.quiver(origin[0], origin[1], 1, 0, axis_scale, 0, **axis_style)
    ax.text(origin[0]-axis_scale*1.3, origin[1], 1, 'D', ha='right', va='center', **label_style)
    ax.text(origin[0]+axis_scale*1.3, origin[1], 1, 'V', ha='left', va='center', **label_style)
    ax.text(origin[0], origin[1]+axis_scale*1.3, 1, 'I', ha='center', va='bottom', **label_style)
    ax.text(0.5, 0, 1, 'Stereo-seq', **label_style)

    # ATAC-layer axes
    draw_axes([-0.48, -0.48], -1, ['UMAP1', 'UMAP2'], axis_scale)
    ax.text(0.5, 0, -1, 'snATAC-seq', **label_style)

    # ================== View optimization ==================
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()

    # Title
    ax.text2D(0.5, 0.9, title, transform=ax.transAxes,
              ha='center', va='top', fontsize=14, fontweight='bold', color='#404040')

    # Legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                       label=cluster, markerfacecolor=color, markersize=10)
                      for cluster, color in target_color_map.items() if cluster in target_subtypes]
    fig.legend(handles=legend_elements, loc='center left', 
               bbox_to_anchor=(0.85, 0.5), title="Cell Types",
               title_fontsize=12, fontsize=10, frameon=False)

    plt.tight_layout()
    plt.close()
    return fig