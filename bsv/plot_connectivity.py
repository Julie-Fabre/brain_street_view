import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def plot_connectivity(experiment_data, allen_atlas_path, output_region,
                      number_of_chunks, number_of_pixels, plane,
                      region_only, smoothing, color_limits, color, normalization_info='unknown',
                      input_regions=None, region_groups=None, experiment_region_info=None,
                      normalize_by_group=False, custom_slices=None, slice_averaging=0,
                      atlas_type='allen', atlas_resolution=10):
    if input_regions is None:
        input_regions = []
    if region_groups is None:
        region_groups = []
    if custom_slices is None:
        custom_slices = []

    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)
    atlas_slice_spacing = atlas_resolution

    # Ensure output_region is a string
    if isinstance(output_region, (list, tuple)):
        output_region = output_region[0]

    # Find structure indices for output region
    curr_plot_structure_idx = find_structure_indices(st, output_region)
    plot_structure_color = get_structure_color(st, curr_plot_structure_idx[0])

    # Handle region groups
    if input_regions:
        n_regions = len(input_regions)
        if region_groups:
            if len(region_groups) != n_regions:
                raise ValueError(f'Length of region_groups ({len(region_groups)}) must match input_regions ({n_regions})')
            unique_groups = sorted(set(region_groups))
            n_region_groups = len(unique_groups)
            region_groups_cell = [[j for j, rg in enumerate(region_groups) if rg == g] for g in unique_groups]
        else:
            n_region_groups = n_regions
            region_groups_cell = [[i] for i in range(n_regions)]
    else:
        n_region_groups = 1
        region_groups_cell = [[0]]

    # Get chunk limits
    # Only use left hemisphere of atlas for structure limits
    half_ml = av.shape[2] // 2
    structure_mask = np.isin(av[:, :, :half_ml], curr_plot_structure_idx)
    ap_vals, _, ml_vals = np.where(structure_mask)

    if custom_slices:
        custom_indices = [s * atlas_slice_spacing for s in custom_slices]
        chunks_region = []
        for ci in custom_indices:
            chunks_region.append(ci - slice_averaging * atlas_slice_spacing)
            chunks_region.append(ci + slice_averaging * atlas_slice_spacing)
        chunks_region = sorted(set(chunks_region))
        number_of_chunks = len(custom_slices)
    else:
        if plane == 'coronal':
            limits = [ap_vals.min(), ap_vals.max()]
        else:
            limits = [ml_vals.min(), ml_vals.max()]
        step = (limits[1] - limits[0]) / number_of_chunks
        chunks_region = [limits[0] + i * step for i in range(number_of_chunks + 1)]

    # Second pass: get ML x DV (coronal) or AP x DV (sagittal) bins
    boundary_projection = [None] * number_of_chunks
    projection_view_bins = [None] * number_of_chunks

    for i_chunk in range(number_of_chunks):
        chunk_start = int(round(chunks_region[i_chunk]))
        chunk_end = int(round(chunks_region[i_chunk + 1]))

        if plane == 'coronal':
            region_area = np.isin(av[chunk_start:chunk_end + 1, :, :half_ml], curr_plot_structure_idx)
            # AP, DV, ML -> get ML, AP, DV coordinates
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
            ap_loc = ap_loc + chunk_start
            x_coords = ml_loc
            y_coords = dv_loc
        else:
            region_area = np.isin(av[:, :, chunk_start:chunk_end + 1], curr_plot_structure_idx)
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
            ml_loc = ml_loc + chunk_start
            x_coords = ap_loc
            y_coords = dv_loc

        if len(x_coords) < 3:
            projection_view_bins[i_chunk] = [np.array([0]), np.array([0])]
            boundary_projection[i_chunk] = np.array([])
            continue

        # Compute boundary using convex hull
        from scipy.spatial import ConvexHull
        points = np.column_stack([x_coords, y_coords])
        try:
            hull = ConvexHull(points)
            hull_indices = np.append(hull.vertices, hull.vertices[0])
        except Exception:
            hull_indices = np.arange(len(x_coords))
        boundary_projection[i_chunk] = hull_indices

        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()

        x_edges = np.linspace(x_min, x_max, number_of_pixels + 1)
        y_edges = np.linspace(y_min, y_max, number_of_pixels + 1)
        projection_view_bins[i_chunk] = [x_edges, y_edges]

    # Collapse hemispheres
    half_slices = experiment_data.shape[2] // 2
    if experiment_data.ndim == 4:
        n_groups = experiment_data.shape[3]
        collapsed = experiment_data[:, :, :half_slices, :] + experiment_data[:, :, ::-1, :][:, :, :half_slices, :]
    else:
        n_groups = 1
        collapsed = experiment_data[:, :, :half_slices] + experiment_data[:, :, ::-1][:, :, :half_slices]

    # Extract projection data for each chunk
    projection_matrix = [None] * number_of_chunks
    for i_chunk in range(number_of_chunks):
        x_edges = projection_view_bins[i_chunk][0]
        y_edges = projection_view_bins[i_chunk][1]

        if 0 in x_edges or 0 in y_edges:
            if n_groups > 1:
                projection_matrix[i_chunk] = np.zeros((number_of_pixels + 1, number_of_pixels + 1, n_groups))
            else:
                projection_matrix[i_chunk] = np.zeros((number_of_pixels + 1, number_of_pixels + 1))
            continue

        this_diff = np.mean(np.diff(chunks_region))

        if plane == 'coronal':
            ap_start = max(0, int(round((chunks_region[i_chunk] - this_diff) / 10)))
            ap_end = min(collapsed.shape[0] - 1, int(round((chunks_region[i_chunk] + this_diff) / 10)))
            y_idx = np.clip((y_edges / 10).astype(int), 0, collapsed.shape[1] - 1)
            x_idx = np.clip((x_edges / 10).astype(int), 0, collapsed.shape[2] - 1 if collapsed.ndim >= 3 else 0)

            if collapsed.ndim == 4:
                data_slice = collapsed[ap_start:ap_end + 1][:, y_idx][:, :, x_idx, :]
                mean_data = np.nanmean(data_slice, axis=0)  # Average over AP
                projtemp = mean_data.transpose(1, 0, 2)  # DV x ML x groups -> ML x DV x groups
            else:
                data_slice = collapsed[ap_start:ap_end + 1][:, y_idx][:, :, x_idx]
                mean_data = np.nanmean(data_slice, axis=0)
                projtemp = mean_data.T
        else:
            ml_start = max(0, int(round((chunks_region[i_chunk] - this_diff) / 10)))
            ml_end = min(collapsed.shape[2] - 1 if collapsed.ndim >= 3 else 0,
                         int(round((chunks_region[i_chunk] + this_diff) / 10)))
            x_idx = np.clip((x_edges / 10).astype(int), 0, collapsed.shape[0] - 1)
            y_idx = np.clip((y_edges / 10).astype(int), 0, collapsed.shape[1] - 1)

            if collapsed.ndim == 4:
                data_slice = collapsed[x_idx][:, y_idx][:, :, ml_start:ml_end + 1, :]
                mean_data = np.nanmean(data_slice, axis=2)
                projtemp = mean_data.transpose(1, 0, 2)
            else:
                data_slice = collapsed[x_idx][:, y_idx][:, :, ml_start:ml_end + 1]
                mean_data = np.nanmean(data_slice, axis=2)
                projtemp = mean_data.T

        projection_matrix[i_chunk] = projtemp

    # Group-wise normalization
    if normalize_by_group and projection_matrix[0].ndim == 3 and projection_matrix[0].shape[2] > 1:
        n_g = projection_matrix[0].shape[2]
        group_max = np.zeros(n_g)
        for i_g in range(n_g):
            for i_chunk in range(number_of_chunks):
                group_max[i_g] = max(group_max[i_g], np.nanmax(projection_matrix[i_chunk][:, :, i_g]))
        for i_chunk in range(number_of_chunks):
            for i_g in range(n_g):
                if group_max[i_g] > 0:
                    projection_matrix[i_chunk][:, :, i_g] /= group_max[i_g]

    # Plot
    fig, axes = plt.subplots(n_region_groups, number_of_chunks,
                              figsize=(3 * number_of_chunks, 3 * n_region_groups),
                              squeeze=False)
    fig.patch.set_facecolor('white')
    fig.canvas.manager.set_window_title('Fluorescence intensity')

    slice_aras = np.zeros(number_of_chunks)

    for i_chunk in range(number_of_chunks):
        # Get boundary for region masking
        chunk_start = int(round(chunks_region[i_chunk]))
        chunk_end = int(round(chunks_region[i_chunk + 1]))

        if plane == 'coronal':
            region_area = np.isin(av[chunk_start:chunk_end + 1, :, :half_ml], curr_plot_structure_idx)
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
        else:
            region_area = np.isin(av[:, :, chunk_start:chunk_end + 1], curr_plot_structure_idx)
            ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))

        x_edges = projection_view_bins[i_chunk][0]
        y_edges = projection_view_bins[i_chunk][1]

        # Build in-polygon mask
        if plane == 'coronal':
            bnd_x, bnd_y = ml_loc, dv_loc
        else:
            bnd_x, bnd_y = ap_loc, dv_loc

        is_in = _build_region_mask(x_edges, y_edges, bnd_x, bnd_y)

        for i_rg in range(n_region_groups):
            ax = axes[i_rg, i_chunk]

            if projection_matrix[i_chunk].ndim == 3 and projection_matrix[i_chunk].shape[2] > i_rg:
                avg_data = projection_matrix[i_chunk][:, :, i_rg]
            else:
                avg_data = projection_matrix[i_chunk] if projection_matrix[i_chunk].ndim == 2 else projection_matrix[i_chunk][:, :, 0]

            # Mask outside region
            masked_data = avg_data.copy()
            masked_data[~is_in] = np.nan

            ax.imshow(masked_data.T, origin='upper' if plane == 'coronal' else 'lower',
                       extent=[x_edges[0], x_edges[-1], y_edges[-1], y_edges[0]] if plane == 'coronal'
                       else [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
                       cmap='gray_r', vmin=0, vmax=1, aspect='equal')
            ax.set_facecolor('0.5')

            # Plot boundary
            if len(bnd_x) >= 3:
                from scipy.spatial import ConvexHull
                pts = np.column_stack([bnd_x, bnd_y])
                try:
                    hull = ConvexHull(pts)
                    hull_pts = pts[np.append(hull.vertices, hull.vertices[0])]
                    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=plot_structure_color, linewidth=2)
                except Exception:
                    pass

            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')

            # ARA slice level
            if custom_slices:
                this_slice_ara = custom_slices[i_chunk]
            else:
                this_slice_ara = int(round(np.nanmean(chunks_region[i_chunk:i_chunk + 2]) / 10))
            slice_aras[i_chunk] = this_slice_ara

            if i_rg == 0:
                if i_chunk == 0:
                    prefix = 'ARA level (cor.): ' if plane == 'coronal' else 'ARA level (sag.): '
                    ax.set_title(f'{prefix}{this_slice_ara}', fontsize=9)
                else:
                    ax.set_title(str(this_slice_ara), fontsize=9)

            if i_chunk == 0 and input_regions:
                regions_in_group = region_groups_cell[i_rg]
                group_names = [input_regions[r] for r in regions_in_group]
                label = '+'.join(group_names)
                ax.text(-0.15, 0.5, label, transform=ax.transAxes,
                        fontweight='bold', fontsize=12, ha='right', va='center', rotation=90)

    plt.tight_layout()
    plt.show(block=False)

    # Build return arrays
    if n_groups <= 1:
        proj_array = np.stack([projection_matrix[s] if projection_matrix[s].ndim == 2
                               else projection_matrix[s][:, :, 0]
                               for s in range(number_of_chunks)], axis=-1)
    else:
        proj_array = np.stack(projection_matrix, axis=-1)

    proj_coords = []
    for s in range(number_of_chunks):
        coords = list(projection_view_bins[s]) + [slice_aras[s] * 10]
        proj_coords.append(coords)

    return proj_array, proj_coords


def _build_region_mask(x_edges, y_edges, bnd_x, bnd_y):
    """Build boolean mask of pixels inside the region boundary."""
    nx = len(x_edges)
    ny = len(y_edges)
    is_in = np.zeros((nx, ny), dtype=bool)

    if len(bnd_x) < 3:
        return is_in

    from scipy.spatial import ConvexHull
    pts = np.column_stack([bnd_x, bnd_y])
    try:
        hull = ConvexHull(pts)
        hull_pts = pts[hull.vertices]
        path = Path(hull_pts)
        for ix in range(nx):
            for iy in range(ny):
                is_in[ix, iy] = path.contains_point([x_edges[ix], y_edges[iy]])
    except Exception:
        is_in[:] = True

    return is_in
