import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from skimage import measure
import pandas as pd
from scipy.spatial import ConvexHull
import brewer2mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch

def load_structure_tree(file_path):
    return pd.read_csv(file_path)

def plot_connectivity(experiment_data, allen_atlas_path, input_region, number_of_chunks, number_of_pixels, plane, region_only, smoothing, color_limits, color):
    # Load Allen atlas
    av = np.load(os.path.join(allen_atlas_path, 'annotation_volume_10um_by_index.npy'))
    st = load_structure_tree(os.path.join(allen_atlas_path, 'structure_tree_safe_2017.csv'))
    atlas_slice_spacing = 10  # 10 um/ slice

    # Current region info
    curr_plot_structure_idx = st[st['acronym'].str.contains(input_region)].index
    keep_struct = [st['acronym'].iloc[i].startswith(input_region) for i in curr_plot_structure_idx]
    curr_plot_structure_idx = curr_plot_structure_idx[keep_struct]

    plot_structure_color = np.array([int(st['color_hex_triplet'].iloc[curr_plot_structure_idx[0]][i:i+2], 16) for i in (0, 2, 4)]) / 255.0

    # Get chunk limits in AP (if coronal) or ML (if sagital)
    structure_limits = np.where(np.isin(av[:, :, :av.shape[2]//2], curr_plot_structure_idx))
    ap_values, _, ml_values = structure_limits

    if plane == 'coronal':
        curr_limits = [np.min(ap_values), np.max(ap_values)]
        chunks_region = np.linspace(curr_limits[0], curr_limits[1], number_of_chunks + 1)
        projection_views = np.tile([0, 1], (number_of_chunks, 1))  # ML x AP
    elif plane == 'sagital':
        curr_limits = [np.min(ml_values), np.max(ml_values)]
        chunks_region = np.linspace(curr_limits[0], curr_limits[1], number_of_chunks + 1)
        projection_views = np.tile([1, 0], (number_of_chunks, 1))  # AP x ML

    # Initialize variables
    boundary_projection = [None] * number_of_chunks
    projection_view_bins = [None] * number_of_chunks
    projection_view_lims = np.full((number_of_chunks, 2, 2), np.nan)

    fig_outline1 = plt.figure('Chunk AP limits' if plane == 'coronal' else 'Chunk ML limits')

    for i_chunk in range(number_of_chunks):
        if plane == 'coronal':
            region_area = np.isin(av[int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1]), :, :av.shape[2]//2], curr_plot_structure_idx).transpose(2, 0, 1)
        else:
            region_area = np.isin(av[:, :, int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1])], curr_plot_structure_idx).transpose(2, 0, 1)

        region_location = np.array(np.where(region_area)).T

        this_chunk_x = region_location[:, projection_views[0, 0]]
        this_chunk_y = region_location[:, projection_views[0, 1]]

        if plane == 'coronal':
            ax = fig_outline1.add_subplot(number_of_chunks, 1, i_chunk + 1)
        else:
            ax = fig_outline1.add_subplot(1, number_of_chunks, i_chunk + 1)

        hull = ConvexHull(np.column_stack((this_chunk_x, this_chunk_y)))
        boundary_projection[i_chunk] = hull.vertices

        if plane == 'coronal':
            ax.plot(region_location[hull.vertices, projection_views[0, 0]], 
                     region_location[hull.vertices, projection_views[0, 1]], 
                     color=plot_structure_color)
        else:
            ax.plot(region_location[hull.vertices, projection_views[0, 1]], 
                     region_location[hull.vertices, projection_views[0, 0]], 
                     color=plot_structure_color)

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        projection_view_lims[i_chunk, 0, :] = ax.get_xlim()
        projection_view_lims[i_chunk, 1, :] = ax.get_ylim()
        projection_view_bins[i_chunk] = [
            np.linspace(projection_view_lims[i_chunk, 0, 0], projection_view_lims[i_chunk, 0, 1], number_of_pixels + 1),
            np.linspace(projection_view_lims[i_chunk, 1, 0], projection_view_lims[i_chunk, 1, 1], number_of_pixels + 1)
        ]

        if plane == 'coronal':
            ax.invert_yaxis()

    prettify_plot(fig_outline1)

    # Get ML x DV values (if coronal) or AP x DV values (if sagital)
    projection_grid_size = [132, 80, 114]
    fig_outline2 = plt.figure('Chunk ML x DV limits' if plane == 'coronal' else 'Chunk AP x DV limits')

    for i_chunk in range(number_of_chunks):
        if plane == 'coronal':
            region_area = np.isin(av[int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1]), :, :av.shape[2]//2], curr_plot_structure_idx).transpose(2, 0, 1)
        else:
            region_area = np.isin(av[:, :, int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1])], curr_plot_structure_idx).transpose(2, 0, 1)

        region_location = np.array(np.where(region_area)).T
        if plane == 'coronal':
            region_location[:, 1] += int(chunks_region[i_chunk])
        else:
            region_location[:, 0] += int(chunks_region[i_chunk])

        this_chunk_x = region_location[:, projection_views[0, 0]]
        this_chunk_dv = region_location[:, projection_views[0, 1]]

        ax = fig_outline2.add_subplot(1, number_of_chunks, i_chunk + 1)
        hull = ConvexHull(np.column_stack((this_chunk_x, this_chunk_dv)))
        boundary_projection[i_chunk] = hull.vertices

        ax.plot(region_location[hull.vertices, projection_views[0, 0]], 
                 region_location[hull.vertices, projection_views[0, 1]], 
                 color=plot_structure_color)

        projection_view_lims[i_chunk, 0, :] = ax.get_xlim()
        projection_view_lims[i_chunk, 1, :] = ax.get_ylim()
        projection_view_bins[i_chunk] = [
            np.linspace(projection_view_lims[i_chunk, 0, 0], projection_view_lims[i_chunk, 0, 1], number_of_pixels + 1),
            np.linspace(projection_view_lims[i_chunk, 1, 0], projection_view_lims[i_chunk, 1, 1], number_of_pixels + 1)
        ]

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    prettify_plot(fig_outline2)

    # Process experiment data
    these_locations_ap_dv_ml = np.zeros(experiment_data.shape[:3] + (experiment_data.shape[3] // 2,))
    for i_group in range(experiment_data.shape[3]):
        these_locations_ap_dv_ml[:, :, :, i_group] = (
            experiment_data[:, :, :projection_grid_size[2]//2, i_group] +
            experiment_data[:, :, projection_grid_size[2]-1:projection_grid_size[2]//2-1:-1, i_group]
        )

    projection_matrix = {}
    for i_chunk in range(number_of_chunks):
        projection_matrix[i_chunk] = np.zeros((number_of_pixels + 1, number_of_pixels + 1, experiment_data.shape[3]))

        x_edges = projection_view_bins[i_chunk][0]
        y_edges = projection_view_bins[i_chunk][1]

        this_diff = np.mean(np.diff(chunks_region))

        if not np.any(x_edges == 0) and not np.any(y_edges == 0):
            if plane == 'coronal':
                proj_temp = np.nanmean(these_locations_ap_dv_ml[
                    int((chunks_region[i_chunk] - this_diff) / 10):int((chunks_region[i_chunk] + this_diff) / 10),
                    np.round(y_edges / 10).astype(int),
                    np.round(x_edges / 10).astype(int),
                    :
                ], axis=0).transpose(1, 0, 2)
            else:
                proj_temp = np.nanmean(these_locations_ap_dv_ml[
                    np.round(x_edges / 10).astype(int),
                    np.round(y_edges / 10).astype(int),
                    int((chunks_region[i_chunk] - this_diff) / 10):int((chunks_region[i_chunk] + this_diff) / 10),
                    :
                ], axis=2).transpose(1, 0, 2)

            projection_matrix[i_chunk] = proj_temp

    # Plot fluorescence intensity
    fig_projection = plt.figure('Fluorescence intensity', figsize=(15, 5 * experiment_data.shape[3]))
    n_groups = experiment_data.shape[3]

    for i_chunk in range(number_of_chunks):
        if plane == 'coronal':
            region_area = np.isin(av[int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1]), :, :av.shape[2]//2], curr_plot_structure_idx).transpose(2, 0, 1)
        else:
            region_area = np.isin(av[:, :, int(chunks_region[i_chunk]):int(chunks_region[i_chunk+1])], curr_plot_structure_idx).transpose(2, 0, 1)

        region_location = np.array(np.where(region_area)).T

        is_in = np.zeros((projection_matrix[i_chunk].shape[0], projection_matrix[i_chunk].shape[1]), dtype=bool)
        path = Path(np.column_stack((region_location[boundary_projection[i_chunk], projection_views[0, 0]], 
                                     region_location[boundary_projection[i_chunk], projection_views[0, 1]])))
        
        for i_pixel_x in range(projection_matrix[i_chunk].shape[0]):
            for i_pixel_y in range(projection_matrix[i_chunk].shape[1]):
                is_in[i_pixel_x, i_pixel_y] = path.contains_point((projection_view_bins[i_chunk][0][i_pixel_x], 
                                                                   projection_view_bins[i_chunk][1][i_pixel_y]))

        for i_group in range(n_groups):
            ax = fig_projection.add_subplot(n_groups, number_of_chunks, i_group * number_of_chunks + i_chunk + 1)

            max_value = np.nanmax([np.nanmax(projection_matrix[c][:, :, i_group]) for c in range(number_of_chunks)])
            this_cmap_limits = [0, max_value]

            binned_array_pixel = projection_matrix[i_chunk][:, :, i_group]
            binned_array_pixel[~is_in] = np.nan

            # Apply smoothing if required
            if smoothing:
                binned_array_pixel_smooth = smooth_data(binned_array_pixel, smoothing)
            else:
                binned_array_pixel_smooth = binned_array_pixel

            im = ax.imshow(binned_array_pixel_smooth.T, extent=[projection_view_bins[i_chunk][0][0], projection_view_bins[i_chunk][0][-1],
                                                                projection_view_bins[i_chunk][1][0], projection_view_bins[i_chunk][1][-1]],
                           cmap='Greys', vmin=this_cmap_limits[0], vmax=this_cmap_limits[1])

            ax.set_facecolor([0.5, 0.5, 0.5])
            im.set_clim(this_cmap_limits)

            ax.plot(region_location[boundary_projection[i_chunk], projection_views[0, 0]], 
                    region_location[boundary_projection[i_chunk], projection_views[0, 1]], 
                    color=plot_structure_color, linewidth=2)

            ax.set_aspect('equal')
            ax.axis('off')

            if i_chunk == 0:
                this_slice_ara = round(np.nanmean(chunks_region[i_chunk:i_chunk+2]) / 10)
                ax.set_title(f"ARA level ({'cor' if plane == 'coronal' else 'sag'}): {this_slice_ara}")
            elif i_chunk == number_of_chunks - 1:
                this_slice_ara = round(np.nanmean(chunks_region[i_chunk:i_chunk+2]) / 10)
                ax.set_title(f"{this_slice_ara}")
            else:
                this_slice_ara = round(np.nanmean(chunks_region[i_chunk:i_chunk+2]) / 10)
                ax.set_title(f"{this_slice_ara}")

    # Set x and y limits
    xlims_region = np.array([ax.get_xlim() for ax in fig_projection.axes[::number_of_chunks]])
    ylims_region = np.array([ax.get_ylim() for ax in fig_projection.axes[::number_of_chunks]])

    diff_xlims_region = np.diff(xlims_region, axis=1).flatten()
    diff_ylims_region = np.diff(ylims_region, axis=1).flatten()

    for i_chunk in range(number_of_chunks):
        for i_group in range(n_groups):
            ax = fig_projection.axes[i_group * number_of_chunks + i_chunk]

            xlims_here = (np.max(diff_xlims_region) - diff_xlims_region[i_chunk]) / 2
            ax.set_xlim([xlims_region[i_chunk, 0] - xlims_here, xlims_region[i_chunk, 1] + xlims_here])

            ylims_here = (np.max(diff_ylims_region) - diff_ylims_region[i_chunk]) / 2
            ax.set_ylim([ylims_region[i_chunk, 0] - ylims_here, ylims_region[i_chunk, 1] + ylims_here])

            shrink_factor_x = diff_xlims_region[i_chunk] / (diff_xlims_region[i_chunk] + xlims_here * 2)

            if i_chunk == 0 and i_group == 0:
                axis_length_mm = 1
                one_pixel_x = np.diff(projection_view_bins[i_chunk][0][1:3])[0]
                one_pixel_x_um = one_pixel_x / 2.5 / shrink_factor_x  # QQ: Why 2.5 here?
                axis_length_atlas_units_x = (axis_length_mm * 1000) / one_pixel_x_um
                add_scale_bar(ax, axis_length_atlas_units_x, 0, f"{axis_length_mm}mm", "", "topLeft", "", "")

    plt.tight_layout()
    return fig_outline1, fig_outline2, fig_projection

def smooth_data(data, smoothing):
    # Implement smoothing here. You might use gaussian_filter from scipy.ndimage, for example.
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(data, sigma=smoothing)

def add_scale_bar(ax, width, height, width_label, height_label, location, color, fontcolor):
    # This is a simplified version. You might need to adjust it based on your specific requirements.
    if location == "topLeft":
        x, y = 0.05, 0.95
    else:
        # Add other locations as needed
        x, y = 0.05, 0.95

    ax.plot([x, x + width], [y, y], color=color or 'k', linewidth=2)
    ax.plot([x, x], [y, y - height], color=color or 'k', linewidth=2)
    
    if width_label:
        ax.text(x + width/2, y, width_label, ha='center', va='bottom', color=fontcolor or 'k')
    if height_label:
        ax.text(x, y - height/2, height_label, ha='right', va='center', color=fontcolor or 'k', rotation=90)

def prettify_plot(fig, options=None):
    # This is a placeholder. Implement according to your specific prettify_plot function.
    for ax in fig.axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

# Main execution
if __name__ == "__main__":
    # Example usage:
    experiment_data = np.random.rand(132, 80, 114, 2)  # Random data for demonstration
    allen_atlas_path = "/path/to/allen/atlas"
    input_region = "VISp"
    number_of_chunks = 5
    number_of_pixels = 100
    plane = "coronal"
    region_only = True
    smoothing = 2
    color_limits = [0, 1]
    color = [1, 0, 0]

    fig_outline1, fig_outline2, fig_projection = plot_connectivity(
        experiment_data, allen_atlas_path, input_region, number_of_chunks, 
        number_of_pixels, plane, region_only, smoothing, color_limits, color
    )

    plt.show()


def plot_connectivity_3d(injection_summary, allen_atlas_path, region_to_plot, color, plot_patch):
    # Load atlas
    av = np.load(os.path.join(allen_atlas_path, 'annotation_volume_10um_by_index.npy'))
    st = pd.read_csv(os.path.join(allen_atlas_path, 'structure_tree_safe_2017.csv'))

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot brain outline/volume
    fig, _ = plot_brain_grid(ax=ax)

    # Plot region
    curr_plot_structure_idx = st[st['acronym'].str.contains(region_to_plot)].index
    region_to_plot_structure_color = np.array([int(st.loc[curr_plot_structure_idx[0], 'color_hex_triplet'][i:i+2], 16) for i in (0, 2, 4)]) / 255.0

    slice_spacing = 10
    structure_3d = measure.marching_cubes(np.transpose(np.isin(av[::slice_spacing, ::slice_spacing, ::slice_spacing], curr_plot_structure_idx), (2, 0, 1)), 0)

    structure_alpha = 0.3
    ax.plot_trisurf(structure_3d[0][:, 0]*slice_spacing, 
                    structure_3d[0][:, 1]*slice_spacing, 
                    structure_3d[0][:, 2]*slice_spacing,
                    triangles=structure_3d[1],
                    color=region_to_plot_structure_color,
                    alpha=structure_alpha)

    # Plot injections
    ax.set_box_aspect((1,1,1))

    max_voxel_x = injection_summary['max_voxel_x'] / 10
    max_voxel_y = injection_summary['max_voxel_y'] / 10
    max_voxel_z = injection_summary['max_voxel_z'] / 10

    # Remove any missing data
    keep_me = max_voxel_x != 0
    injection_structure_idx = [st.index[st['id'] == id].tolist()[0] for id in injection_summary['structure_id']]

    plot_all_colors = False  # leave false for now
    if plot_all_colors:
        hex_triplets = st.loc[injection_structure_idx[keep_me], 'color_hex_triplet']
        rgb_matrix = np.array([[int(color[i:i+2], 16) for i in (0, 2, 4)] for color in hex_triplets]) / 255.0
    else:
        rgb_matrix = np.array([int(st.loc[injection_structure_idx[0], 'color_hex_triplet'][i:i+2], 16) for i in (0, 2, 4)]) / 255.0

    dot_size = 10**3
    scatter = ax.scatter(max_voxel_x[keep_me], 
                         max_voxel_z[keep_me], 
                         max_voxel_y[keep_me],
                         s=injection_summary['projection_volume'][keep_me] * dot_size,
                         c=rgb_matrix,
                         alpha=0.5)

    plt.gcf().set_facecolor('w')

    # Uncomment to save as rotating video
    # from celluloid import Camera
    # camera = Camera(fig)
    # for angle in range(0, 360, 10):
    #     ax.view_init(elev=10, azim=angle)
    #     camera.snap()
    # animation = camera.animate()
    # animation.save('rotating_brain.mp4')

    plt.show()

def plot_brain_grid(brain_grid_data=None, ax=None, brain_figure=None, black_brain=False):
    """
    Plot the wire mesh data loaded from brainGridData.npy.

    Parameters:
    brain_grid_data (numpy.ndarray): The brain grid data to plot. If None, it will be loaded from file.
    ax (matplotlib.axes.Axes): The axes to plot on. If None, a new figure and axes will be created.
    brain_figure (matplotlib.figure.Figure): The figure to plot on. If None, a new figure will be created.
    black_brain (bool): If True, use a black background for the brain plot.

    Returns:
    tuple: (figure, plot_handle)
    """
    if brain_grid_data is None:
        # Get the directory of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        brain_grid_data = np.load(os.path.join(current_dir, 'brainGridData.npy'))

    bp = brain_grid_data.astype(float)
    # When saved to uint16, NaN's become zeros. 
    # There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were
    bp[np.sum(bp, axis=1) == 0, :] = np.nan

    if ax is None:
        if brain_figure is None:
            brain_figure = plt.figure(figsize=(10, 10), name='Brain View')
        ax = brain_figure.add_subplot(111, projection='3d')

    if not black_brain:
        h = ax.plot(bp[:, 0], bp[:, 1], bp[:, 2], color=(0, 0, 0, 0.3))
    else:
        h = ax.plot(bp[:, 0], bp[:, 1], bp[:, 2], color=(0.7, 0.7, 0.7, 0.3))
        ax.get_figure().patch.set_facecolor('k')

    ax.set_zlim(ax.get_zlim()[::-1])  # Reverse the Z-axis
    ax.set_box_aspect((1, 1, 1))  # Equal aspect ratio
    ax.set_axis_off()

    plt.tight_layout()

    return ax.get_figure(), h

# Example usage:
# fig, plot_handle = plot_brain_grid()
# plt.show()