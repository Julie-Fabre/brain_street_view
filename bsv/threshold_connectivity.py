import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.ndimage import label, gaussian_filter
from matplotlib.path import Path

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def threshold_connectivity(experiment_data, allen_atlas_path, input_region,
                           number_of_chunks, number_of_pixels, plane,
                           region_only, smoothing, color_limits, color,
                           threshold, threshold_method='absolute',
                           normalization_method='none', data_fetch_normalization='unknown',
                           atlas_type='allen', atlas_resolution=10):
    """Threshold projection density and visualize significant signals.

    Parameters
    ----------
    experiment_data : numpy.ndarray
        Projection density array from :func:`fetch_connectivity_data`.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    input_region : str
        Target region acronym (e.g. ``'CP'``).
    number_of_chunks : int
        Number of slices.
    number_of_pixels : int
        Pixel resolution per slice panel.
    plane : str
        ``'coronal'`` or ``'sagittal'``.
    region_only : bool
        Mask to the target region boundary.
    smoothing : float
        Gaussian smoothing sigma in pixels.
    color_limits : str or list
        Colour scale specification.
    color : list or None
        RGB colour(s).
    threshold : float
        Threshold value (interpretation depends on *threshold_method*).
    threshold_method : str, optional
        ``'absolute'``, ``'percentile'``, ``'zscore'``, or ``'relative'``.
    normalization_method : str, optional
        Additional normalization: ``'none'``, ``'region'``, ``'zscore'``,
        or ``'robust'``.
    data_fetch_normalization : str, optional
        Label for the normalization applied during data fetch.
    atlas_type : str, optional
        Atlas type (default ``'allen'``).
    atlas_resolution : int, optional
        Atlas resolution in micrometres (10 or 20).

    Returns
    -------
    proj_array : numpy.ndarray
        Thresholded binary projection matrix.
    proj_coords : list
        Coordinate information for each slice.
    """
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)
    atlas_slice_spacing = atlas_resolution

    if isinstance(input_region, (list, tuple)):
        input_region = input_region[0]

    curr_idx = find_structure_indices(st, input_region)
    plot_structure_color = get_structure_color(st, curr_idx[0])

    # Get chunk limits
    half_ml = av.shape[2] // 2
    structure_mask = np.isin(av[:, :, :half_ml], curr_idx)
    ap_vals, _, ml_vals = np.where(structure_mask)

    if plane == 'coronal':
        limits = [ap_vals.min(), ap_vals.max()]
    else:
        limits = [ml_vals.min(), ml_vals.max()]
    step = (limits[1] - limits[0]) / number_of_chunks
    chunks_region = [limits[0] + i * step for i in range(number_of_chunks + 1)]

    if plane == 'coronal':
        projection_views = 'ml_dv'
    else:
        projection_views = 'ap_dv'

    # Build bin edges for each chunk
    boundary_projection = [None] * number_of_chunks
    projection_view_bins = [None] * number_of_chunks

    for i_chunk in range(number_of_chunks):
        chunk_start = int(round(chunks_region[i_chunk]))
        chunk_end = int(round(chunks_region[i_chunk + 1]))

        if plane == 'coronal':
            region_area = np.isin(av[chunk_start:chunk_end + 1, :, :half_ml], curr_idx)
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
            ap_loc += chunk_start
            x_coords, y_coords = ml_loc, dv_loc
        else:
            region_area = np.isin(av[:, :, chunk_start:chunk_end + 1], curr_idx)
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
            ml_loc += chunk_start
            x_coords, y_coords = ap_loc, dv_loc

        if len(x_coords) < 3:
            projection_view_bins[i_chunk] = [np.array([0]), np.array([0])]
            boundary_projection[i_chunk] = (x_coords, y_coords)
            continue

        boundary_projection[i_chunk] = (x_coords, y_coords)
        x_edges = np.linspace(x_coords.min(), x_coords.max(), number_of_pixels + 1)
        y_edges = np.linspace(y_coords.min(), y_coords.max(), number_of_pixels + 1)
        projection_view_bins[i_chunk] = [x_edges, y_edges]

    # Collapse hemispheres
    n_groups = experiment_data.shape[3] if experiment_data.ndim == 4 else 1
    half_slices = experiment_data.shape[2] // 2
    if experiment_data.ndim == 4:
        collapsed = np.zeros((*experiment_data.shape[:2], half_slices, n_groups))
        for g in range(n_groups):
            collapsed[:, :, :, g] = (experiment_data[:, :, :half_slices, g] +
                                     experiment_data[:, :, -1:half_slices - 1:-1, g])
    else:
        collapsed = experiment_data[:, :, :half_slices] + experiment_data[:, :, -1:half_slices - 1:-1]

    # Conversion factor: atlas voxels per projection grid voxel
    atlas_to_grid = 100 / atlas_resolution

    # Extract projection data
    projection_matrix = [None] * number_of_chunks
    for i_chunk in range(number_of_chunks):
        x_edges = projection_view_bins[i_chunk][0]
        y_edges = projection_view_bins[i_chunk][1]

        if 0 in x_edges or 0 in y_edges:
            projection_matrix[i_chunk] = np.zeros((number_of_pixels + 1, number_of_pixels + 1))
            continue

        this_diff = np.mean(np.diff(chunks_region))

        if plane == 'coronal':
            ap_s = max(0, int(round((chunks_region[i_chunk] - this_diff) / atlas_to_grid)))
            ap_e = min(collapsed.shape[0] - 1, int(round((chunks_region[i_chunk] + this_diff) / atlas_to_grid)))
            y_idx = np.clip((y_edges / atlas_to_grid).astype(int), 0, collapsed.shape[1] - 1)
            x_idx = np.clip((x_edges / atlas_to_grid).astype(int), 0, collapsed.shape[2] - 1 if collapsed.ndim >= 3 else 0)
            if collapsed.ndim == 4:
                data_slice = collapsed[ap_s:ap_e + 1][:, y_idx][:, :, x_idx, :]
                mean_data = np.nanmean(data_slice, axis=0)
                projtemp = mean_data.transpose(1, 0, 2)
            else:
                data_slice = collapsed[ap_s:ap_e + 1][:, y_idx][:, :, x_idx]
                mean_data = np.nanmean(data_slice, axis=0)
                projtemp = mean_data.T
        else:
            ml_s = max(0, int(round((chunks_region[i_chunk] - this_diff) / atlas_to_grid)))
            ml_e = min(collapsed.shape[2] - 1 if collapsed.ndim >= 3 else 0,
                       int(round((chunks_region[i_chunk] + this_diff) / atlas_to_grid)))
            x_idx = np.clip((x_edges / atlas_to_grid).astype(int), 0, collapsed.shape[0] - 1)
            y_idx = np.clip((y_edges / atlas_to_grid).astype(int), 0, collapsed.shape[1] - 1)
            if collapsed.ndim == 4:
                data_slice = collapsed[x_idx][:, y_idx][:, :, ml_s:ml_e + 1, :]
                mean_data = np.nanmean(data_slice, axis=2)
                projtemp = mean_data.transpose(1, 0, 2)
            else:
                data_slice = collapsed[x_idx][:, y_idx][:, :, ml_s:ml_e + 1]
                mean_data = np.nanmean(data_slice, axis=2)
                projtemp = mean_data.T

        projection_matrix[i_chunk] = projtemp

    # Step 1: Global 0-1 normalization
    all_data = np.concatenate([pm.ravel() for pm in projection_matrix])
    all_data = all_data[~np.isnan(all_data)]
    global_min, global_max = all_data.min(), all_data.max()
    print(f'Global data range: {global_min:.6f} to {global_max:.6f}')

    if global_max > global_min:
        for i_chunk in range(number_of_chunks):
            projection_matrix[i_chunk] = (projection_matrix[i_chunk] - global_min) / (global_max - global_min)

    # Step 2: Additional normalization
    if normalization_method != 'none':
        for i_chunk in range(number_of_chunks):
            pm = projection_matrix[i_chunk]
            if pm.ndim == 3:
                for g in range(pm.shape[2]):
                    projection_matrix[i_chunk][:, :, g] = _apply_normalization(pm[:, :, g], normalization_method)
            else:
                projection_matrix[i_chunk] = _apply_normalization(pm, normalization_method)

    # Calculate adaptive threshold
    all_data = np.concatenate([pm.ravel() for pm in projection_matrix])
    all_data = all_data[~np.isnan(all_data)]

    if threshold_method == 'percentile':
        adaptive_threshold = np.percentile(all_data, threshold)
        print(f'Calculated {threshold}th percentile threshold: {adaptive_threshold:.4f}')
    elif threshold_method == 'zscore':
        adaptive_threshold = np.mean(all_data) + threshold * np.std(all_data)
        print(f'Calculated {threshold}-sigma threshold: {adaptive_threshold:.4f}')
    elif threshold_method == 'relative':
        adaptive_threshold = threshold * np.max(all_data)
        print(f'Calculated relative threshold: {adaptive_threshold:.4f}')
    else:
        adaptive_threshold = threshold
        print(f'Using absolute threshold: {adaptive_threshold:.4f}')

    # Statistics
    n_above = np.sum(all_data > adaptive_threshold)
    print(f'Voxels above threshold: {n_above} ({100 * n_above / len(all_data):.1f}%)')
    print(f'Mean +/- SD: {np.mean(all_data):.4f} +/- {np.std(all_data):.4f}')

    # Store original for second figure
    original_matrix = [pm.copy() for pm in projection_matrix]

    # Apply threshold (binary)
    for i_chunk in range(number_of_chunks):
        projection_matrix[i_chunk] = (projection_matrix[i_chunk] > adaptive_threshold).astype(float)

    # Plot thresholded figure
    fig_thresh, axes_t = plt.subplots(max(1, n_groups), number_of_chunks,
                                       figsize=(3 * number_of_chunks, 3 * max(1, n_groups)),
                                       squeeze=False)
    fig_thresh.patch.set_facecolor('white')
    fig_thresh.canvas.manager.set_window_title('Thresholded Connectivity')

    # Plot original + boundary figure
    fig_orig, axes_o = plt.subplots(max(1, n_groups), number_of_chunks,
                                     figsize=(3 * number_of_chunks, 3 * max(1, n_groups)),
                                     squeeze=False)
    fig_orig.patch.set_facecolor('white')
    fig_orig.canvas.manager.set_window_title('Original with Threshold Boundary')

    slice_aras = np.zeros(number_of_chunks)

    for i_chunk in range(number_of_chunks):
        bnd_x, bnd_y = boundary_projection[i_chunk]
        x_edges = projection_view_bins[i_chunk][0]
        y_edges = projection_view_bins[i_chunk][1]

        is_in = _build_mask(x_edges, y_edges, bnd_x, bnd_y)

        this_slice_ara = int(round(np.nanmean(chunks_region[i_chunk:i_chunk + 2]) / atlas_to_grid))
        slice_aras[i_chunk] = this_slice_ara

        for i_group in range(max(1, n_groups)):
            # Get data
            if projection_matrix[i_chunk].ndim == 3:
                thresh_data = projection_matrix[i_chunk][:, :, i_group].copy()
                orig_data = original_matrix[i_chunk][:, :, i_group].copy()
            else:
                thresh_data = projection_matrix[i_chunk].copy()
                orig_data = original_matrix[i_chunk].copy()

            thresh_data[~is_in] = np.nan
            orig_data[~is_in] = np.nan

            # Apply Gaussian smoothing if requested
            if smoothing and smoothing > 0:
                for name in ['thresh', 'orig']:
                    arr = thresh_data if name == 'thresh' else orig_data
                    nan_mask = np.isnan(arr)
                    temp = gaussian_filter(np.nan_to_num(arr, nan=0.0), sigma=smoothing)
                    temp[nan_mask] = np.nan
                    if name == 'thresh':
                        thresh_data = temp
                    else:
                        orig_data = temp

            # Thresholded plot
            ax_t = axes_t[i_group, i_chunk]
            ax_t.imshow(thresh_data.T, origin='upper' if plane == 'coronal' else 'lower',
                         extent=[x_edges[0], x_edges[-1], y_edges[-1], y_edges[0]],
                         cmap='gray_r', vmin=0, vmax=1, aspect='equal')
            ax_t.set_facecolor('0.5')
            _plot_boundary(ax_t, bnd_x, bnd_y, plot_structure_color)
            ax_t.axis('off')
            if i_group == 0:
                ax_t.set_title(f'ARA {this_slice_ara}', fontsize=9)

            # Original + boundary plot
            ax_o = axes_o[i_group, i_chunk]
            ax_o.imshow(orig_data.T, origin='upper' if plane == 'coronal' else 'lower',
                         extent=[x_edges[0], x_edges[-1], y_edges[-1], y_edges[0]],
                         cmap='gray_r', vmin=0, vmax=1, aspect='equal')
            ax_o.set_facecolor('0.5')
            _plot_boundary(ax_o, bnd_x, bnd_y, plot_structure_color)

            # Plot threshold boundary on original
            if projection_matrix[i_chunk].ndim == 3:
                thresh_mask = original_matrix[i_chunk][:, :, 0] > adaptive_threshold
            else:
                thresh_mask = original_matrix[i_chunk] > adaptive_threshold
            _plot_threshold_boundary(ax_o, thresh_mask, x_edges, y_edges)

            ax_o.axis('off')
            if i_group == 0:
                prefix = 'Original - ARA ' if i_chunk == 0 else ''
                ax_o.set_title(f'{prefix}{this_slice_ara}', fontsize=9)

    for fig in [fig_thresh, fig_orig]:
        fig.tight_layout()

    plt.show(block=False)

    # Build return arrays
    if n_groups <= 1:
        proj_array = np.stack([pm if pm.ndim == 2 else pm[:, :, 0]
                               for pm in projection_matrix], axis=-1)
    else:
        proj_array = np.stack(projection_matrix, axis=-1)

    proj_coords = []
    for s in range(number_of_chunks):
        coords = list(projection_view_bins[s]) + [slice_aras[s] * 10]
        proj_coords.append(coords)

    return proj_array, proj_coords


def _apply_normalization(data, method):
    if method == 'region':
        mn, mx = np.nanmin(data), np.nanmax(data)
        return (data - mn) / (mx - mn) if mx > mn else data
    elif method == 'zscore':
        m, s = np.nanmean(data), np.nanstd(data)
        return (data - m) / s if s > 0 else data
    elif method == 'robust':
        med = np.nanmedian(data)
        q75, q25 = np.nanpercentile(data, 75), np.nanpercentile(data, 25)
        iqr = q75 - q25
        return (data - med) / iqr if iqr > 0 else data
    return data


def _build_mask(x_edges, y_edges, bnd_x, bnd_y):
    nx, ny = len(x_edges), len(y_edges)
    is_in = np.zeros((nx, ny), dtype=bool)
    if len(bnd_x) < 3:
        return is_in
    try:
        hull = ConvexHull(np.column_stack([bnd_x, bnd_y]))
        path = Path(np.column_stack([bnd_x, bnd_y])[hull.vertices])
        for ix in range(nx):
            for iy in range(ny):
                is_in[ix, iy] = path.contains_point([x_edges[ix], y_edges[iy]])
    except Exception:
        is_in[:] = True
    return is_in


def _plot_boundary(ax, bnd_x, bnd_y, color):
    if len(bnd_x) < 3:
        return
    try:
        hull = ConvexHull(np.column_stack([bnd_x, bnd_y]))
        pts = np.column_stack([bnd_x, bnd_y])[np.append(hull.vertices, hull.vertices[0])]
        ax.plot(pts[:, 0], pts[:, 1], color=color, linewidth=2)
    except Exception:
        pass


def _plot_threshold_boundary(ax, mask, x_edges, y_edges):
    if not mask.any():
        return
    from scipy.ndimage import binary_dilation
    boundary = binary_dilation(mask) ^ mask
    by, bx = np.where(boundary.T)
    if len(bx) < 2:
        return
    # Convert to coordinate space
    if len(x_edges) > 1 and len(y_edges) > 1:
        x_scale = (x_edges[-1] - x_edges[0]) / (len(x_edges) - 1)
        y_scale = (y_edges[-1] - y_edges[0]) / (len(y_edges) - 1)
        ax.plot(x_edges[0] + bx * x_scale, y_edges[0] + by * y_scale,
                'r.', markersize=2)
