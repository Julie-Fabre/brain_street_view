import os
import json
import numpy as np
import pandas as pd

from .fetch_connectivity_summary import fetch_connectivity_summary, load_injection_summary
from .fetch_connectivity_images import fetch_connectivity_images, PROJECTION_GRID_SIZE
from .atlas_utils import load_projection_info


def fetch_connectivity_data(experiment_ids, save_location, file_name,
                            normalization_method, subtract_other_hemisphere,
                            grouping_method='', allen_atlas_path='',
                            load_all=False, input_regions=None, region_groups=None,
                            export_metadata=True, reload=False,
                            atlas_type='allen', atlas_resolution=10):
    """Download and cache projection density maps for a set of experiments.

    Parameters
    ----------
    experiment_ids : list of int
        Experiment IDs returned by :func:`find_connectivity_experiments`.
    save_location : str
        Local directory for cached downloads.
    file_name : str
        Base name for the cached metadata CSV. ``''`` to skip caching.
    normalization_method : str
        ``'none'`` or ``'injectionIntensity'`` (divide by injection volume).
    subtract_other_hemisphere : bool
        Subtract the contralateral hemisphere signal.
    grouping_method : str, optional
        Group experiments by ``'AP'``, ``'ML'``, ``'DV'``, or ``''`` (no
        grouping).
    allen_atlas_path : str, optional
        Path to the Allen CCF atlas directory.
    load_all : bool, optional
        If True, also return individual (non-averaged) projection volumes.
    input_regions : list of str, optional
        Region acronyms for region-based grouping.
    region_groups : list of int, optional
        Group assignment per region (same length as *input_regions*).
    export_metadata : bool, optional
        Write a per-experiment metadata CSV.
    reload : bool, optional
        Re-download even if cached files exist.
    atlas_type : str, optional
        Atlas type (default ``'allen'``).
    atlas_resolution : int, optional
        Atlas resolution in micrometres (10 or 20).

    Returns
    -------
    combined_projection : numpy.ndarray
        Averaged projection density array (AP x DV x ML x groups).
    combined_injection_info : dict
        Aggregated injection metadata across experiments.
    individual_projections : numpy.ndarray or None
        Per-experiment volumes when *load_all* is True, else None.
    experiment_region_info : dict
        Per-experiment region metadata.
    """
    if not grouping_method:
        grouping_method = 'NaN'
    if input_regions is None:
        input_regions = []
    if region_groups is None:
        region_groups = []

    os.makedirs(save_location, exist_ok=True)

    # Load projection info CSV
    projection_info = load_projection_info()

    n_exp = len(experiment_ids)
    primary_structure_AP = np.zeros(n_exp)
    primary_structure_DV = np.zeros(n_exp)
    primary_structure_ML = np.zeros(n_exp)
    primary_structure_ID = np.zeros(n_exp, dtype=int)
    primary_structure_abbreviation = [''] * n_exp

    # Collect all injection info
    combined_injection_info = {}

    print(f'Loading {n_exp} experiments...')

    for i, exp_id in enumerate(experiment_ids):
        save_dir = os.path.join(save_location, str(exp_id))
        os.makedirs(save_dir, exist_ok=True)

        # Fetch summary if needed
        summary_path = os.path.join(save_dir, 'injectionSummary_all.json')
        if not os.path.exists(summary_path):
            status = fetch_connectivity_summary(exp_id, save_dir)
            if not status:
                continue

        injection_info = load_injection_summary(save_dir)

        # Accumulate injection info
        for entry in injection_info:
            for key, val in entry.items():
                if key not in combined_injection_info:
                    combined_injection_info[key] = []
                combined_injection_info[key].append(val)

        # Get metadata from projection info CSV
        row = projection_info[projection_info['id'] == exp_id]
        if len(row) > 0:
            row = row.iloc[0]
            primary_structure_ID[i] = row['structure_id']
            primary_structure_abbreviation[i] = row['structure_abbrev']
            coords_str = str(row['injection_coordinates'])
            coords = [float(x) for x in coords_str.strip('[]').split(',')]
            if len(coords) >= 3:
                primary_structure_AP[i] = coords[0]
                primary_structure_DV[i] = coords[1]
                primary_structure_ML[i] = coords[2]

    # Print experiment summary
    exp_info = projection_info[projection_info['id'].isin(experiment_ids)]
    print(f'\nEXPERIMENT SUMMARY')
    print(f'Total experiments loaded: {n_exp}')

    print(f'\nMouse genotype distribution:')
    genotypes = exp_info['transgenic_line'].fillna('Wild-type').replace({'': 'Wild-type', '""': 'Wild-type'})
    for gt, count in genotypes.value_counts().items():
        print(f'  {gt}: {count} experiments')

    print(f'\nBrain region distribution:')
    for reg, count in exp_info['structure_abbrev'].value_counts().items():
        print(f'  {reg}: {count} experiments')

    # Grouping
    if grouping_method in ('AP', 'ML', 'DV'):
        coord_map = {'AP': primary_structure_AP, 'ML': primary_structure_ML, 'DV': primary_structure_DV}
        vals = coord_map[grouping_method]
        sorted_idx = np.argsort(vals)
        unique_vals = np.unique(vals[sorted_idx])
        group_id = np.zeros(n_exp, dtype=int)
        for j, v in enumerate(unique_vals):
            group_id[vals == v] = j + 1
        n_groups = len(unique_vals)
    elif grouping_method == 'NaN' or not grouping_method:
        group_id = np.ones(n_exp, dtype=int)
        n_groups = 1
    else:
        print(f'Warning: grouping method "{grouping_method}" not recognized - skipping grouping.')
        group_id = np.ones(n_exp, dtype=int)
        n_groups = 1

    # Region-based grouping override
    if input_regions and region_groups:
        unique_rg = sorted(set(region_groups))
        n_region_groups = len(unique_rg)
        rg_cell = {g: [j for j, rg in enumerate(region_groups) if rg == g] for g in unique_rg}

        region_based_group_id = np.zeros(n_exp, dtype=int)
        for i in range(n_exp):
            abbrev = primary_structure_abbreviation[i]
            region_idx = -1
            for j, reg in enumerate(input_regions):
                if abbrev == reg or abbrev.startswith(reg):
                    region_idx = j
                    break
            if region_idx >= 0:
                for g_idx, g in enumerate(unique_rg):
                    if region_idx in rg_cell[g]:
                        region_based_group_id[i] = g_idx + 1
                        break

        n_groups = n_region_groups
        group_id = region_based_group_id
        print(f'Region-based grouping: {np.sum(group_id > 0)} experiments mapped to {n_groups} groups')

    # Initialize combined projection
    combined_projection = np.zeros((*PROJECTION_GRID_SIZE, n_groups))
    individual_projections = np.zeros((*PROJECTION_GRID_SIZE, n_exp)) if load_all else None

    # Load raw images and accumulate
    print('Getting raw images...')
    for i, exp_id in enumerate(experiment_ids):
        curr_group = group_id[i]
        if curr_group == 0:
            continue

        raw_path = os.path.join(save_location, str(exp_id), 'density.raw')
        if not os.path.exists(raw_path):
            status = fetch_connectivity_images(exp_id, os.path.join(save_location, str(exp_id)))
            if not status:
                continue

        # Read density.raw - MATLAB reads column-major (Fortran order)
        experiment_projection = np.fromfile(raw_path, dtype='<f4').reshape(PROJECTION_GRID_SIZE, order='F')

        # Normalize by injection volume (mm^3).
        # The density.raw values are already fractional (0-1), so dividing
        # by injection volume makes experiments with different injection
        # sizes comparable.
        if normalization_method in ('injectionVolume', 'injectionIntensity'):
            exp_indices = np.array(combined_injection_info.get('experimentID', [])) == exp_id
            struct_indices = exp_indices & (np.array(combined_injection_info.get('structure_id', [])) == 997)

            hem_ids = np.array(combined_injection_info.get('hemisphere_id', []))
            hem3 = struct_indices & (hem_ids == 3)

            vol_arr = np.array(combined_injection_info.get('projection_volume', []), dtype=float)
            inj_vol = vol_arr[hem3].max() if hem3.any() else 0

            if inj_vol == 0:
                print(f'Warning: Invalid injection volume for experiment {exp_id}, skipping normalization')
                inj_vol = 1

            experiment_projection = experiment_projection / inj_vol

        # Subtract other hemisphere
        if subtract_other_hemisphere:
            injection_info = load_injection_summary(os.path.join(save_location, str(exp_id)))
            max_voxel_z = injection_info[0].get('max_voxel_z', 57)
            tmp = np.zeros_like(experiment_projection)
            if max_voxel_z <= 57:  # left hemisphere
                for ml in range(57):
                    tmp[:, :, ml] = experiment_projection[:, :, ml] - experiment_projection[:, :, 113 - ml]
            else:  # right hemisphere
                for ml in range(57):
                    tmp[:, :, 113 - ml] = experiment_projection[:, :, 113 - ml] - experiment_projection[:, :, ml]
            experiment_projection = tmp

        if load_all:
            individual_projections[:, :, :, i] = experiment_projection
        combined_projection[:, :, :, curr_group - 1] += experiment_projection

    # Average by group count
    for g in range(n_groups):
        count = np.sum(group_id == g + 1)
        if count > 0:
            combined_projection[:, :, :, g] /= count

    # Build experiment region info
    experiment_region_info = {
        'experimentIDs': experiment_ids,
        'primaryStructure_ID': primary_structure_ID,
        'primaryStructure_abbreviation': primary_structure_abbreviation,
        'primaryStructure_AP': primary_structure_AP,
        'primaryStructure_DV': primary_structure_DV,
        'primaryStructure_ML': primary_structure_ML,
    }

    # Export metadata CSV
    if export_metadata:
        csv_base = file_name if file_name else 'connectivity_data'
        meta = pd.DataFrame()
        meta['ExperimentID'] = experiment_ids

        exp_info_ordered = projection_info[projection_info['id'].isin(experiment_ids)]
        idx_map = {eid: i for i, eid in enumerate(experiment_ids)}
        exp_info_ordered = exp_info_ordered.copy()
        exp_info_ordered['_sort'] = exp_info_ordered['id'].map(idx_map)
        exp_info_ordered = exp_info_ordered.sort_values('_sort').drop(columns='_sort')

        meta['MouseLine'] = exp_info_ordered['transgenic_line'].fillna('Wild-type').replace(
            {'': 'Wild-type', '""': 'Wild-type'}).values
        meta['InjectionRegion'] = primary_structure_abbreviation
        meta['InjectionRegionID'] = primary_structure_ID
        meta['InjectionAP'] = primary_structure_AP
        meta['InjectionDV'] = primary_structure_DV
        meta['InjectionML'] = primary_structure_ML

        # Injection volume/intensity from combined info
        inj_volumes = np.zeros(n_exp)
        inj_intensities = np.zeros(n_exp)
        exp_id_arr = np.array(combined_injection_info.get('experimentID', []))
        struct_id_arr = np.array(combined_injection_info.get('structure_id', []))
        hem_id_arr = np.array(combined_injection_info.get('hemisphere_id', []))
        pv_arr = np.array(combined_injection_info.get('projection_volume', []), dtype=float)
        spi_arr = np.array(combined_injection_info.get('sum_pixel_intensity', []), dtype=float)

        for j, eid in enumerate(experiment_ids):
            mask = (exp_id_arr == eid) & (struct_id_arr == 997) & (hem_id_arr == 3)
            if mask.any():
                inj_volumes[j] = pv_arr[mask].max()
                inj_intensities[j] = spi_arr[mask].max()

        meta['InjectionVolume'] = inj_volumes
        meta['InjectionIntensity'] = inj_intensities

        if 'sex' in exp_info_ordered.columns:
            meta['Sex'] = exp_info_ordered['sex'].values
        if 'age' in exp_info_ordered.columns:
            meta['Age'] = exp_info_ordered['age'].values

        csv_path = os.path.join(save_location, f'{csv_base}_metadata.csv')
        meta.to_csv(csv_path, index=False)
        print(f'Metadata exported to: {csv_path}')

    print('Data processing complete.')
    return combined_projection, combined_injection_info, individual_projections, experiment_region_info
