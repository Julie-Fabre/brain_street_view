import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import bsv


def test_analyze_cp_subregions(atlas_path, save_location, test_experiment_ids):
    """Test CP subregion analysis with real data."""
    v2_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF_v2'
    if not os.path.exists(v2_atlas_path):
        import pytest
        pytest.skip('v2 atlas not found')

    combined, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False, allen_atlas_path=atlas_path)

    # Get projection data first
    proj_array, proj_coords = bsv.plot_connectivity(
        combined, atlas_path, 'CP', 5, 10, 'coronal',
        True, 2, 'global', None, 'none')

    sub_results, global_results = bsv.analyze_cp_subregions(
        proj_array, proj_coords, v2_atlas_path)

    assert 'subregion_ids' in sub_results
    assert 'mean_intensities' in sub_results
    assert 'subregion_ids' in global_results
    assert 'mean_intensities' in global_results

    # Should find some subregions with data
    n_with_data = np.sum(~np.isnan(global_results['mean_intensities']))
    print(f'Subregions with data: {n_with_data}')
    assert n_with_data > 0


def test_combined_analysis(atlas_path, save_location, test_experiment_ids):
    """Test the combined connectivity + subregion analysis."""
    v2_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF_v2'
    if not os.path.exists(v2_atlas_path):
        import pytest
        pytest.skip('v2 atlas not found')

    combined, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False, allen_atlas_path=atlas_path)

    proj_array, proj_coords, sub_results, global_results = bsv.plot_connectivity_with_subregion_analysis(
        combined, atlas_path, v2_atlas_path, 'CP',
        5, 10, 'coronal', True, 2, 'global', None, 'none')

    assert proj_array is not None
    assert sub_results is not None
    assert global_results is not None
    print(f'Projection shape: {proj_array.shape}')
