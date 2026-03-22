import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def plot_connectivity_3d(injection_summary, allen_atlas_path, region_to_plot,
                         color=None, plot_patch=True,
                         atlas_type='allen', atlas_resolution=10):
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)

    # Find region and get color
    curr_idx = find_structure_indices(st, region_to_plot)
    region_color = get_structure_color(st, curr_idx[0])
    slice_spacing = atlas_resolution

    # Subsample and create isosurface
    av_sub = av[::slice_spacing, ::slice_spacing, ::slice_spacing]
    region_mask = np.isin(av_sub, curr_idx).transpose(2, 0, 1)  # ML, AP, DV

    from skimage.measure import marching_cubes
    try:
        verts, faces, _, _ = marching_cubes(region_mask.astype(float), level=0.5)
        verts = verts * slice_spacing
    except Exception:
        print('Warning: Could not generate isosurface')
        verts, faces = np.zeros((0, 3)), np.zeros((0, 3), dtype=int)

    fig = plt.figure(figsize=(10, 8))
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, projection='3d')

    # Plot region surface
    if len(faces) > 0:
        mesh = Poly3DCollection(verts[faces], alpha=0.3, facecolor=region_color, edgecolor='none')
        ax.add_collection3d(mesh)

    # Plot injection sites
    max_voxel_x = np.array(injection_summary.get('max_voxel_x', [])) / atlas_resolution
    max_voxel_y = np.array(injection_summary.get('max_voxel_y', [])) / atlas_resolution
    max_voxel_z = np.array(injection_summary.get('max_voxel_z', [])) / atlas_resolution
    proj_vol = np.array(injection_summary.get('projection_volume', []))

    keep = max_voxel_x != 0
    if keep.any():
        # Get injection color from first structure
        struct_ids = np.array(injection_summary.get('structure_id', []))
        first_struct = struct_ids[keep][0]
        matching = st.index[st['id'] == first_struct].tolist()
        if matching:
            rgb = get_structure_color(st, matching[0] + 1)  # convert 0-based pandas idx to 1-based
        else:
            rgb = region_color

        dot_size = proj_vol[keep] * 1e3
        dot_size = np.clip(dot_size, 5, 500)

        ax.scatter(max_voxel_x[keep], max_voxel_z[keep], max_voxel_y[keep],
                   s=dot_size, c=[rgb], alpha=0.5, edgecolors=[rgb])

    ax.set_xlabel('AP')
    ax.set_ylabel('ML')
    ax.set_zlabel('DV')
    ax.set_box_aspect([1, 1, 1])

    plt.tight_layout()
    plt.show(block=False)
