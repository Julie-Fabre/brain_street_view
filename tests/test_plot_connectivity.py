import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for testing
import bsv


def test_plot_connectivity_basic(atlas_path, save_location, test_experiment_ids):
    """Test basic 2D plotting with pre-fetched data."""
    # Load data without normalization
    combined, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False, allen_atlas_path=atlas_path)

    proj_array, proj_coords = bsv.plot_connectivity(
        combined, atlas_path, 'CP', 5, 10, 'coronal',
        True, 2, 'global', None, 'none')

    assert proj_array.ndim >= 2
    assert proj_array.shape[-1] == 5  # number_of_chunks
    assert len(proj_coords) == 5
    # Each coord entry should have [x_edges, y_edges, ara_slice]
    for c in proj_coords:
        assert len(c) == 3
    print(f'Projection array shape: {proj_array.shape}')
    print(f'Projection max: {proj_array.max():.6f}')


def test_plot_connectivity_sagittal(atlas_path, save_location, test_experiment_ids):
    """Test sagittal plane plotting."""
    combined, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False, allen_atlas_path=atlas_path)

    proj_array, proj_coords = bsv.plot_connectivity(
        combined, atlas_path, 'CP', 3, 8, 'sagital',
        True, 2, 'global', None, 'none')

    assert proj_array.ndim >= 2
    assert proj_array.shape[-1] == 3
    print(f'Sagittal projection shape: {proj_array.shape}')


def test_threshold_connectivity(atlas_path, save_location, test_experiment_ids):
    """Test thresholding with different methods."""
    combined, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False, allen_atlas_path=atlas_path)

    # Percentile threshold
    proj_array, proj_coords = bsv.threshold_connectivity(
        combined, atlas_path, 'CP', 5, 10, 'coronal',
        True, 2, 'global', None,
        threshold=90, threshold_method='percentile',
        normalization_method='none')

    assert proj_array.ndim >= 2
    # Thresholded data should be binary (0 or 1), except NaN
    vals = proj_array[~np.isnan(proj_array)]
    unique_vals = np.unique(vals)
    assert all(v in [0.0, 1.0] for v in unique_vals)
    print(f'Threshold array shape: {proj_array.shape}')
    print(f'Unique values: {unique_vals}')
