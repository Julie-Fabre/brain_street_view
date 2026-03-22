import numpy as np
import bsv
from bsv.atlas_utils import load_atlas, find_structure_indices, get_structure_color, load_projection_info


def test_load_atlas(atlas_path):
    av, st = load_atlas(atlas_path)
    # Allen 10um atlas should be [1320, 800, 1140]
    assert av.ndim == 3
    assert av.shape[0] > 1000  # AP
    assert av.shape[1] > 500   # DV
    assert av.shape[2] > 1000  # ML
    assert len(st) > 500  # many structures


def test_find_structure_indices(atlas_path):
    _, st = load_atlas(atlas_path)
    # CP (Caudoputamen) should exist
    idx = find_structure_indices(st, 'CP')
    assert len(idx) > 0
    # Returned indices are 1-based; check all returned acronyms start with 'CP'
    for i in idx:
        assert st.loc[i - 1, 'acronym'].startswith('CP')


def test_get_structure_color(atlas_path):
    _, st = load_atlas(atlas_path)
    idx = find_structure_indices(st, 'CP')
    color = get_structure_color(st, idx[0])
    assert color.shape == (3,)
    assert np.all(color >= 0) and np.all(color <= 1)


def test_load_projection_info():
    info = load_projection_info()
    assert len(info) > 100
    assert 'id' in info.columns
    assert 'structure_abbrev' in info.columns  # hyphen->underscore from CSV
    assert 'injection_coordinates' in info.columns
