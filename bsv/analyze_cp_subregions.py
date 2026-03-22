import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def analyze_cp_subregions(projection_matrix_array, projection_matrix_coordinates_ara,
                           allen_atlas_path_v2, output_slices=None, input_regions=None,
                           region_groups=None, save_csv_path=''):
    if input_regions is None:
        input_regions = []
    if region_groups is None:
        region_groups = []

    print('\n=== CP AND NAc SUBREGION ANALYSIS ===')

    # Load v2 atlas
    av_v2 = np.load(os.path.join(allen_atlas_path_v2, 'annotation_volume_v2_20um_by_index.npy'))
    ontology_v2 = pd.read_csv(os.path.join(allen_atlas_path_v2, 'UnifiedAtlas_Label_ontology_v2.csv'))

    print(f'Atlas dimensions: {av_v2.shape}')
    print(f'Ontology entries: {len(ontology_v2)}')

    # Find CP and NAc subregions
    cp_mask = (ontology_v2['name'].str.contains('Caudoputamen', case=False, na=False) |
               ontology_v2['acronym'].str.contains('CP', case=False, na=False))
    nac_mask = (ontology_v2['name'].str.contains('Nucleus accumbens', case=False, na=False) |
                ontology_v2['name'].str.contains('accumbens', case=False, na=False) |
                ontology_v2['acronym'].str.contains('ACB', case=False, na=False) |
                ontology_v2['acronym'].str.contains('NAc', case=False, na=False))
    striatum_mask = cp_mask | nac_mask
    striatum = ontology_v2[striatum_mask].reset_index(drop=True)
    n_subregions = len(striatum)
    print(f'Found {n_subregions} striatum subregions (CP + NAc)')

    cp_sub = ontology_v2[cp_mask]
    nac_sub = ontology_v2[nac_mask]
    print(f'\nCP subregions ({len(cp_sub)}):')
    for _, row in cp_sub.iterrows():
        print(f'  {row["id"]}: {row["name"]} ({row["acronym"]})')
    print(f'\nNAc subregions ({len(nac_sub)}):')
    for _, row in nac_sub.iterrows():
        print(f'  {row["id"]}: {row["name"]} ({row["acronym"]})')

    # Determine data structure
    data = projection_matrix_array
    n_coords = len(projection_matrix_coordinates_ara)

    if data.ndim == 4:
        n_actual_slices = data.shape[3]
        n_groups = data.shape[2]
    elif data.ndim == 3:
        third_dim = data.shape[2]
        if n_coords == 1 and third_dim > 1:
            n_actual_slices = 1
            n_groups = third_dim
        elif n_coords == third_dim:
            n_actual_slices = third_dim
            n_groups = 1
        elif n_coords < third_dim and n_coords > 1:
            n_actual_slices = 1
            n_groups = third_dim
        elif third_dim > 1:
            n_actual_slices = 1
            n_groups = third_dim
        else:
            n_actual_slices = third_dim
            n_groups = 1
    else:
        raise ValueError(f'Unexpected data dimensions: {data.shape}')

    print(f'nActualSlices={n_actual_slices}, nGroups={n_groups}')

    if output_slices is None:
        output_slices = list(range(n_actual_slices))
    else:
        output_slices = [s for s in output_slices if s < n_actual_slices]
    n_slices = len(output_slices)

    use_shared_coords = (len(projection_matrix_coordinates_ara) == 1 and n_actual_slices > 1)
    resolution_factor = 2  # 20um / 10um

    # Initialize results
    mean_intensities = np.full((n_subregions, n_groups * n_slices), np.nan)
    voxel_counts = np.zeros((n_subregions, n_groups * n_slices), dtype=int)
    global_intensity_sums = np.zeros(n_subregions)
    global_voxel_counts = np.zeros(n_subregions, dtype=int)

    # Process slices
    for i_slice, slice_idx in enumerate(output_slices):
        if use_shared_coords:
            coords = projection_matrix_coordinates_ara[0]
        elif slice_idx < len(projection_matrix_coordinates_ara):
            coords = projection_matrix_coordinates_ara[slice_idx]
        else:
            continue

        if not coords or len(coords) < 3 or coords[2] is None:
            continue

        ara_coordinate = coords[2]
        atlas_slice_idx = int(round(ara_coordinate / 2))

        if atlas_slice_idx <= 0 or atlas_slice_idx >= av_v2.shape[0]:
            continue

        atlas_slice = av_v2[atlas_slice_idx, :, :]
        x_coords = np.array(coords[0])
        y_coords = np.array(coords[1])
        x_coords_20um = x_coords / resolution_factor
        y_coords_20um = y_coords / resolution_factor

        for i_group in range(n_groups):
            # Get slice data
            if n_groups > 1 and n_actual_slices == 1:
                slice_data = data[:, :, i_group]
            elif data.ndim == 4:
                slice_data = data[:, :, i_group, slice_idx]
            else:
                slice_data = data[:, :, slice_idx]

            for i_sub in range(n_subregions):
                subregion_id = striatum.iloc[i_sub]['id']
                subregion_mask = (atlas_slice == subregion_id)
                if not subregion_mask.any():
                    continue

                atlas_y, atlas_x = np.where(subregion_mask)
                intensities = []

                for v in range(len(atlas_x)):
                    ax, ay = atlas_x[v], atlas_y[v]
                    if (ax < x_coords_20um.min() or ax > x_coords_20um.max() or
                            ay < y_coords_20um.min() or ay > y_coords_20um.max()):
                        continue

                    dx = int(round(np.clip(np.interp(ax, x_coords_20um, np.arange(len(x_coords))),
                                           0, slice_data.shape[0] - 1)))
                    dy = int(round(np.clip(np.interp(ay, y_coords_20um, np.arange(len(y_coords))),
                                           0, slice_data.shape[1] - 1)))

                    val = slice_data[dx, dy]
                    if not np.isnan(val) and val > 0:
                        intensities.append(val)

                if intensities:
                    if n_groups > 1 and n_actual_slices == 1:
                        result_idx = i_group
                    else:
                        result_idx = i_group * n_slices + i_slice
                    if result_idx < mean_intensities.shape[1]:
                        mean_intensities[i_sub, result_idx] = np.mean(intensities)
                        voxel_counts[i_sub, result_idx] = len(intensities)
                    global_intensity_sums[i_sub] += sum(intensities)
                    global_voxel_counts[i_sub] += len(intensities)

    # Global averages
    global_means = np.full(n_subregions, np.nan)
    for i_sub in range(n_subregions):
        if global_voxel_counts[i_sub] > 0:
            global_means[i_sub] = global_intensity_sums[i_sub] / global_voxel_counts[i_sub]

    # Build results
    subregion_results = {
        'slice_numbers': output_slices,
        'subregion_ids': striatum['id'].values,
        'subregion_names': striatum['name'].tolist(),
        'subregion_acronyms': striatum['acronym'].tolist(),
        'nGroups': n_groups,
        'mean_intensities': mean_intensities,
        'voxel_counts': voxel_counts,
        'coordinates_ARA': projection_matrix_coordinates_ara,
    }
    global_results = {
        'subregion_ids': striatum['id'].values,
        'subregion_names': striatum['name'].tolist(),
        'subregion_acronyms': striatum['acronym'].tolist(),
        'mean_intensities': global_means,
        'total_voxel_counts': global_voxel_counts,
    }

    # Summary
    print(f'\n=== SUMMARY ===')
    print(f'Slices analyzed: {n_slices}')
    print(f'Groups analyzed: {n_groups}')
    with_data = np.sum(np.any(~np.isnan(mean_intensities), axis=1))
    print(f'Subregions with data: {with_data}')

    if with_data > 0:
        sort_idx = np.argsort(-np.nan_to_num(global_means, nan=-np.inf))
        print('\nTop 5 subregions by global mean intensity:')
        for rank, idx in enumerate(sort_idx[:5]):
            if not np.isnan(global_means[idx]):
                print(f'  {rank + 1}. {striatum.iloc[idx]["acronym"]}: {global_means[idx]:.4f}')

    # Bar plots
    has_data_idx = np.where(np.any(~np.isnan(mean_intensities), axis=1))[0]
    if len(has_data_idx) > 0:
        n_cols = min(3, n_groups)
        n_rows = int(np.ceil(n_groups / n_cols))
        fig, axes_arr = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows), squeeze=False)
        fig.suptitle('Mean Fluorescence Intensity per Striatum Subregion by Group',
                     fontsize=14, fontweight='bold')

        all_vals = mean_intensities[has_data_idx].ravel()
        all_vals = all_vals[~np.isnan(all_vals)]
        y_max = all_vals.max() * 1.1 if len(all_vals) > 0 else 1

        for i_group in range(n_groups):
            ax = axes_arr[i_group // n_cols, i_group % n_cols]

            if n_groups > 1 and n_actual_slices == 1:
                group_data = mean_intensities[has_data_idx, i_group]
            else:
                cols = slice(i_group * n_slices, (i_group + 1) * n_slices)
                gi = mean_intensities[has_data_idx][:, cols]
                gv = voxel_counts[has_data_idx][:, cols]
                total_i = np.nansum(gi * gv, axis=1)
                total_v = np.nansum(gv, axis=1)
                group_data = np.where(total_v > 0, total_i / total_v, np.nan)

            valid = ~np.isnan(group_data)
            names = [striatum.iloc[idx]['acronym'] for idx in has_data_idx[valid]]
            vals = group_data[valid]

            if len(vals) > 0:
                ax.bar(range(len(vals)), vals, color=[0.2, 0.6, 0.8])
                ax.set_xticks(range(len(names)))
                ax.set_xticklabels(names, rotation=45, ha='right')
                ax.set_ylim(0, y_max)
                ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

            # Title
            if input_regions and region_groups:
                unique_rg = sorted(set(region_groups))
                if i_group < len(unique_rg):
                    g = unique_rg[i_group]
                    reg_indices = [j for j, rg in enumerate(region_groups) if rg == g]
                    group_title = '+'.join([input_regions[j] for j in reg_indices])
                else:
                    group_title = f'Group {i_group + 1}'
            else:
                group_title = f'Group {i_group + 1}'

            ax.set_title(group_title, fontweight='bold')
            ax.set_ylabel('Mean Fluorescence Intensity')
            ax.set_xlabel('Striatum Subregions')

        # Hide unused subplots
        for idx in range(n_groups, n_rows * n_cols):
            axes_arr[idx // n_cols, idx % n_cols].set_visible(False)

        fig.set_facecolor('white')
        plt.tight_layout()
        plt.show(block=False)

    # CSV export
    if save_csv_path:
        summary = pd.DataFrame({
            'SubregionID': striatum['id'].values,
            'SubregionName': striatum['name'].values,
            'SubregionAcronym': striatum['acronym'].values,
            'GlobalMeanIntensity': global_means,
            'TotalVoxelCount': global_voxel_counts,
        })
        summary = summary.sort_values('GlobalMeanIntensity', ascending=False, na_position='last')
        summary.to_csv(save_csv_path, index=False)
        print(f'Results exported to: {save_csv_path}')

        # Non-zero version
        nonzero = summary[summary['GlobalMeanIntensity'].notna() & (summary['GlobalMeanIntensity'] > 0)]
        if len(nonzero) > 0:
            base, ext = os.path.splitext(save_csv_path)
            nonzero.to_csv(f'{base}_nonzero{ext}', index=False)

    return subregion_results, global_results
