import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def plot_connectivity_3d(injection_summary, allen_atlas_path, region_to_plot,
                         color=None, plot_patch=True,
                         atlas_type='allen', atlas_resolution=10):
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)

    curr_idx = find_structure_indices(st, region_to_plot)
    region_color = get_structure_color(st, curr_idx[0])
    slice_spacing = atlas_resolution

    # Subsample atlas: av is [AP, DV, ML]
    av_sub = av[::slice_spacing, ::slice_spacing, ::slice_spacing]

    # MATLAB: permute(ismember(av_sub, idx), [3,1,2]) -> [ML, AP, DV]
    # Then isosurface returns vertices as (col=AP, row=ML, page=DV)
    # In Python: marching_cubes returns (dim0, dim1, dim2)
    # So we feed [ML, AP, DV] and swap cols 0,1 to get (AP, ML, DV)
    region_mask = np.isin(av_sub, curr_idx).transpose(2, 0, 1).astype(float)  # [ML, AP, DV]

    from skimage.measure import marching_cubes
    try:
        verts, faces, _, _ = marching_cubes(region_mask, level=0.5)
        # verts are (ML, AP, DV) — swap to (AP, ML, DV) to match MATLAB
        verts = verts[:, [1, 0, 2]] * slice_spacing
    except Exception:
        print('Warning: Could not generate isosurface')
        verts, faces = np.zeros((0, 3)), np.zeros((0, 3), dtype=int)

    fig = plt.figure(figsize=(10, 8))
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, projection='3d')

    # Plot region surface — vertices are now (AP, ML, DV)
    if len(faces) > 0:
        mesh = Poly3DCollection(verts[faces], alpha=0.3, facecolor=region_color, edgecolor='none')
        ax.add_collection3d(mesh)

    # Injection coordinates: x=AP, y=DV, z=ML (in um)
    # Divide by atlas_resolution to get voxel units matching the isosurface
    max_voxel_x = np.array(injection_summary.get('max_voxel_x', []), dtype=float) / atlas_resolution
    max_voxel_y = np.array(injection_summary.get('max_voxel_y', []), dtype=float) / atlas_resolution
    max_voxel_z = np.array(injection_summary.get('max_voxel_z', []), dtype=float) / atlas_resolution
    proj_vol = np.array(injection_summary.get('projection_volume', []), dtype=float)

    keep = max_voxel_x != 0
    if keep.any():
        struct_ids = np.array(injection_summary.get('structure_id', []))
        first_struct = struct_ids[keep][0]
        matching = st.index[st['id'] == first_struct].tolist()
        if matching:
            rgb = get_structure_color(st, matching[0] + 1)
        else:
            rgb = region_color

        dot_size = proj_vol[keep] * 1e3
        dot_size = np.clip(dot_size, 5, 500)

        # MATLAB: scatter3(x=AP, y=ML, z=DV) = scatter3(max_voxel_x, max_voxel_z, max_voxel_y)
        ax.scatter(max_voxel_x[keep], max_voxel_z[keep], max_voxel_y[keep],
                   s=dot_size, c=[rgb], alpha=0.5, edgecolors=[rgb])

    # Set axis limits from isosurface to ensure both are visible
    if len(verts) > 0:
        for i, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
            margin = (verts[:, i].max() - verts[:, i].min()) * 0.1
            setter(verts[:, i].min() - margin, verts[:, i].max() + margin)

    # Invert z-axis so dorsal (small DV values) is at top
    ax.invert_zaxis()

    ax.set_xlabel('AP (anterior → posterior)')
    ax.set_ylabel('ML (right → left)')
    ax.set_zlabel('DV (dorsal → ventral)')
    ax.set_box_aspect([1, 1, 1])

    plt.tight_layout()
    plt.show(block=False)
