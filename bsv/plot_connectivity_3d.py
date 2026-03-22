import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
try:
    from IPython.display import HTML
except ImportError:
    HTML = None

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def plot_connectivity_3d(injection_summary, allen_atlas_path, region_to_plot,
                         color=None, plot_patch=True, animate=True,
                         atlas_type='allen', atlas_resolution=10):
    """Render injection sites and a target region as a 3D isosurface.

    Parameters
    ----------
    injection_summary : dict
        Injection metadata from :func:`fetch_connectivity_data`.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    region_to_plot : str
        Target region acronym (e.g. ``'CP'``).
    color : list, optional
        RGB colour(s) for injection dots.
    plot_patch : bool, optional
        Render the region as a solid isosurface (True) or grid (False).
    animate : bool, optional
        Create a rotating animation.
    atlas_type : str, optional
        Atlas type (default ``'allen'``).
    atlas_resolution : int, optional
        Atlas resolution in micrometres (10 or 20).

    Returns
    -------
    matplotlib.animation.FuncAnimation or matplotlib.figure.Figure
        Animation object if *animate* is True, otherwise the figure.
    """
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)

    curr_idx = find_structure_indices(st, region_to_plot)
    region_color = get_structure_color(st, curr_idx[0])
    slice_spacing = atlas_resolution

    # Subsample atlas: av is [AP, DV, ML]
    av_sub = av[::slice_spacing, ::slice_spacing, ::slice_spacing]

    # MATLAB: permute([AP,DV,ML], [3,1,2]) -> [ML, AP, DV]
    # isosurface returns (col=AP, row=ML, page=DV)
    # marching_cubes returns (dim0, dim1, dim2), so swap 0,1 to get (AP, ML, DV)
    region_mask = np.isin(av_sub, curr_idx).transpose(2, 0, 1).astype(float)

    from skimage.measure import marching_cubes
    try:
        verts, faces, _, _ = marching_cubes(region_mask, level=0.5)
        verts = verts[:, [1, 0, 2]] * slice_spacing  # -> (AP, ML, DV)
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

    # Injection coordinates (um), divide by atlas_resolution to match isosurface
    max_voxel_x = np.array(injection_summary.get('max_voxel_x', []), dtype=float) / atlas_resolution
    max_voxel_y = np.array(injection_summary.get('max_voxel_y', []), dtype=float) / atlas_resolution
    max_voxel_z = np.array(injection_summary.get('max_voxel_z', []), dtype=float) / atlas_resolution
    proj_vol = np.array(injection_summary.get('projection_volume', []), dtype=float)

    keep = max_voxel_x != 0
    if keep.any():
        struct_ids = np.array(injection_summary.get('structure_id', []))
        first_struct = struct_ids[keep][0]
        matching = st.index[st['id'] == first_struct].tolist()
        rgb = get_structure_color(st, matching[0] + 1) if matching else region_color

        dot_size = proj_vol[keep] * 1e3
        dot_size = np.clip(dot_size, 5, 500)

        # MATLAB: scatter3(AP, ML, DV)
        ax.scatter(max_voxel_x[keep], max_voxel_z[keep], max_voxel_y[keep],
                   s=dot_size, c=[rgb], alpha=0.5, edgecolors=[rgb])

    # Set axis limits from isosurface
    if len(verts) > 0:
        for i, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
            margin = (verts[:, i].max() - verts[:, i].min()) * 0.1
            setter(verts[:, i].min() - margin, verts[:, i].max() + margin)

    # Dorsal at top
    ax.invert_zaxis()

    # Clean look: no grid, no axes
    ax.grid(False)
    ax.set_axis_off()
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    if animate:
        def _rotate(frame):
            ax.view_init(elev=10, azim=frame)
            return []

        anim = FuncAnimation(fig, _rotate, frames=np.arange(0, 360, 2),
                              interval=50, blit=False)
        plt.close(fig)
        if HTML is not None:
            return HTML(anim.to_jshtml())
        return anim
    else:
        plt.tight_layout()
        plt.show(block=False)
        return fig
