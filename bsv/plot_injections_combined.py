import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.path import Path

from .atlas_utils import load_atlas, find_structure_indices, get_structure_color


def plot_injections_combined(experiment_imgs, allen_atlas_path, input_regions,
                              number_of_slices, number_of_pixels, plane,
                              region_only, smoothing, color_limits, color,
                              normalization_method, experiment_region_info=None,
                              atlas_type='allen', atlas_resolution=10):
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)
    n_regions = len(input_regions)

    fig, axes = plt.subplots(n_regions, number_of_slices,
                              figsize=(3 * number_of_slices, 3 * n_regions),
                              squeeze=False)
    fig.patch.set_facecolor('white')
    fig.suptitle('Injection Sites by Region', fontsize=14, fontweight='bold')

    # Collapse hemispheres
    half = experiment_imgs.shape[2] // 2
    collapsed = experiment_imgs[:, :, :half] + experiment_imgs[:, :, ::-1][:, :, :half]
    global_vmax = np.nanmax(collapsed) if np.nanmax(collapsed) > 0 else 1

    for i_reg in range(n_regions):
        region = input_regions[i_reg]
        curr_idx = find_structure_indices(st, region)
        if not curr_idx:
            print(f'Warning: No structure found for region {region}')
            continue
        region_color = get_structure_color(st, curr_idx[0])

        # Get structure limits
        half_ml = av.shape[2] // 2
        structure_mask = np.isin(av[:, :, :half_ml], curr_idx)
        ap_vals, _, ml_vals = np.where(structure_mask)

        if plane == 'coronal':
            limits = [ap_vals.min(), ap_vals.max()]
        else:
            limits = [ml_vals.min(), ml_vals.max()]
        step = (limits[1] - limits[0]) / number_of_slices
        chunks = [limits[0] + j * step for j in range(number_of_slices + 1)]

        for i_chunk in range(number_of_slices):
            ax = axes[i_reg, i_chunk]
            chunk_start = int(round(chunks[i_chunk]))
            chunk_end = int(round(chunks[i_chunk + 1]))

            if plane == 'coronal':
                region_area = np.isin(av[chunk_start:chunk_end + 1, :, :half_ml], curr_idx)
                ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
                ap_loc += chunk_start
                bnd_x, bnd_y = ml_loc, dv_loc
            else:
                region_area = np.isin(av[:, :, chunk_start:chunk_end + 1], curr_idx)
                ml_loc, ap_loc, dv_loc = np.where(region_area.transpose(2, 0, 1))
                ml_loc += chunk_start
                bnd_x, bnd_y = ap_loc, dv_loc

            if len(bnd_x) < 3:
                ax.axis('off')
                continue

            x_range = [bnd_x.min(), bnd_x.max()]
            y_range = [bnd_y.min(), bnd_y.max()]
            if x_range[1] - x_range[0] == 0 or y_range[1] - y_range[0] == 0:
                ax.axis('off')
                continue

            x_edges = np.linspace(x_range[0], x_range[1], number_of_pixels + 1)
            y_edges = np.linspace(y_range[0], y_range[1], number_of_pixels + 1)

            this_diff = np.mean(np.diff(chunks))
            if plane == 'coronal':
                ap_s = max(0, int(round((chunks[i_chunk] - this_diff) / 10)))
                ap_e = min(collapsed.shape[0] - 1, int(round((chunks[i_chunk] + this_diff) / 10)))
                y_idx = np.clip((y_edges / 10).astype(int), 0, collapsed.shape[1] - 1)
                x_idx = np.clip((x_edges / 10).astype(int), 0, collapsed.shape[2] - 1)
                data_slice = collapsed[ap_s:ap_e + 1][:, y_idx][:, :, x_idx]
                mean_data = np.nanmean(data_slice, axis=0)
                projtemp = mean_data.T
            else:
                ax.axis('off')
                continue

            # Mask outside region
            try:
                hull = ConvexHull(np.column_stack([bnd_x, bnd_y]))
                hull_pts = np.column_stack([bnd_x, bnd_y])[hull.vertices]
                path = Path(hull_pts)
                mask = np.zeros(projtemp.shape[:2], dtype=bool)
                for ix in range(len(x_edges)):
                    for iy in range(len(y_edges)):
                        mask[ix, iy] = path.contains_point([x_edges[ix], y_edges[iy]])
                projtemp[~mask] = np.nan
            except Exception:
                pass

            ax.imshow(projtemp.T, origin='upper' if plane == 'coronal' else 'lower',
                       extent=[x_edges[0], x_edges[-1], y_edges[-1], y_edges[0]],
                       cmap='gray_r', vmin=0, vmax=global_vmax, aspect='equal')
            ax.set_facecolor('0.5')

            # Plot boundary
            try:
                hull = ConvexHull(np.column_stack([bnd_x, bnd_y]))
                hull_pts_plot = np.column_stack([bnd_x, bnd_y])[np.append(hull.vertices, hull.vertices[0])]
                ax.plot(hull_pts_plot[:, 0], hull_pts_plot[:, 1], color=region_color, linewidth=2)
            except Exception:
                pass

            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')

            if i_chunk == 0:
                ax.set_ylabel(region, fontweight='bold', fontsize=10, rotation=90)
                ax.axis('on')
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)

            if i_reg == 0:
                slice_ara = int(round(np.mean(chunks[i_chunk:i_chunk + 2]) / 10))
                ax.set_title(f'ARA {slice_ara}', fontsize=9)

    plt.tight_layout()
    plt.show(block=False)
