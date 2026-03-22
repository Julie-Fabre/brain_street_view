"""Tests that run without atlas files (CI-friendly)."""


def test_import():
    import bsv
    assert hasattr(bsv, 'find_connectivity_experiments')
    assert hasattr(bsv, 'fetch_connectivity_data')
    assert hasattr(bsv, 'plot_connectivity')
    assert hasattr(bsv, 'plot_connectivity_3d')
    assert hasattr(bsv, 'threshold_connectivity')
    assert hasattr(bsv, 'analyze_cp_subregions')


def test_projection_info():
    from bsv.atlas_utils import load_projection_info
    info = load_projection_info()
    assert len(info) > 100
    assert 'id' in info.columns
    assert 'structure_abbrev' in info.columns
    assert 'injection_coordinates' in info.columns


def test_find_experiments_api():
    """Smoke test against live Allen API (single small query)."""
    import bsv
    ids = bsv.find_connectivity_experiments(['VISp'])
    assert len(ids) > 0
    assert all(isinstance(x, int) for x in ids)
