import os
import json
import numpy as np
import bsv
from bsv.fetch_connectivity_summary import fetch_connectivity_summary, load_injection_summary
from bsv.fetch_connectivity_images import fetch_connectivity_images, PROJECTION_GRID_SIZE


def test_fetch_summary(save_location, test_experiment_ids):
    exp_id = test_experiment_ids[0]
    save_dir = os.path.join(save_location, str(exp_id))

    # May already be cached, or may need to fetch
    json_path = os.path.join(save_dir, 'injectionSummary_all.json')
    if not os.path.exists(json_path):
        status = fetch_connectivity_summary(exp_id, save_dir)
        if not status:
            import pytest
            pytest.skip('API unreachable, cannot test fetch_summary')

    data = load_injection_summary(save_dir)
    assert len(data) > 0
    # Check expected fields
    entry = data[0]
    assert 'structure_id' in entry
    assert 'hemisphere_id' in entry
    assert 'projection_volume' in entry
    assert 'max_voxel_x' in entry
    assert entry['experimentID'] == exp_id
    print(f'Summary has {len(data)} injection entries')


def test_fetch_images(save_location, test_experiment_ids):
    exp_id = test_experiment_ids[0]
    save_dir = os.path.join(save_location, str(exp_id))
    status = fetch_connectivity_images(exp_id, save_dir)
    assert status

    raw_path = os.path.join(save_dir, 'density.raw')
    assert os.path.exists(raw_path)

    # Read and verify dimensions
    data = np.fromfile(raw_path, dtype='<f4')
    assert data.size == np.prod(PROJECTION_GRID_SIZE)
    data = data.reshape(PROJECTION_GRID_SIZE, order='F')
    assert data.shape == PROJECTION_GRID_SIZE
    assert data.max() > 0  # should have some signal
    print(f'Density data range: {data.min():.6f} to {data.max():.6f}')


def test_fetch_connectivity_data(save_location, atlas_path, test_experiment_ids):
    # Test with no normalization first (injectionIntensity divides by huge sum_pixel_intensity)
    combined, inj_info, _, region_info = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'none', False,
        allen_atlas_path=atlas_path, load_all=False)

    assert combined.shape[:3] == PROJECTION_GRID_SIZE
    assert combined.shape[3] >= 1
    assert combined.max() > 0
    assert 'experimentIDs' in region_info
    assert len(region_info['experimentIDs']) == len(test_experiment_ids)
    print(f'Combined projection shape: {combined.shape}')
    print(f'Max projection value: {combined.max():.6f}')

    # Also test with injectionIntensity - should still run without error
    combined_norm, _, _, _ = bsv.fetch_connectivity_data(
        test_experiment_ids, save_location, '',
        'injectionIntensity', False,
        allen_atlas_path=atlas_path, load_all=False)
    assert combined_norm.shape[:3] == PROJECTION_GRID_SIZE
    print(f'Normalized max: {combined_norm.max():.10f}')
