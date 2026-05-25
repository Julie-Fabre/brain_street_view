import numpy as np
import pandas as pd
import pytest
import bsv
from bsv.atlas_utils import (
    load_atlas,
    get_atlas_files,
    find_structure_indices,
    get_structure_color,
    load_projection_info,
)


# ---------------------------------------------------------------------------
# Pure-logic tests: no atlas files or network required (run everywhere / in CI)
# ---------------------------------------------------------------------------

def _toy_structure_tree():
    """A minimal structure tree with the columns the helpers depend on."""
    return pd.DataFrame({
        'acronym': ['root', 'CP', 'CPu', 'CTX', 'VISp'],
        'color_hex_triplet': ['FFFFFF', '98D6F9', '00FF00', 'ABCDEF', '0080FF'],
    })


def test_get_atlas_files_allen_10um():
    annotation, structure = get_atlas_files('allen', 10)
    assert annotation == 'annotation_volume_10um_by_index.npy'
    assert structure == 'structure_tree_safe_2017.csv'


def test_get_atlas_files_allen_20um():
    annotation, structure = get_atlas_files('allen', 20)
    assert annotation.endswith('.npy')
    assert structure.endswith('.csv')


def test_get_atlas_files_custom_atlas():
    annotation, structure = get_atlas_files('myatlas', 25)
    assert annotation == 'myatlas_annotation_25um.npy'
    assert structure == 'myatlas_structure_tree.csv'


def test_get_atlas_files_bad_resolution():
    with pytest.raises(ValueError):
        get_atlas_files('allen', 99)


def test_find_structure_indices_prefix_match():
    st = _toy_structure_tree()
    # 'CP' should match 'CP' (row 1) and 'CPu' (row 2); returned 1-based.
    assert find_structure_indices(st, 'CP') == [2, 3]


def test_find_structure_indices_exact_prefix():
    st = _toy_structure_tree()
    assert find_structure_indices(st, 'VISp') == [5]


def test_find_structure_indices_no_match():
    st = _toy_structure_tree()
    assert find_structure_indices(st, 'ZZZ') == []


def test_get_structure_color_hex_to_rgb():
    st = _toy_structure_tree()
    # 1-based index 2 -> row 1 -> 'CP' -> '98D6F9'
    color = get_structure_color(st, 2)
    assert color.shape == (3,)
    assert np.allclose(color, [0x98 / 255, 0xD6 / 255, 0xF9 / 255])


def test_get_structure_color_pads_short_hex():
    st = pd.DataFrame({'acronym': ['x'], 'color_hex_triplet': ['FF']})
    # 'FF' -> zero-padded to '0000FF' -> pure blue
    assert np.allclose(get_structure_color(st, 1), [0.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# Atlas-dependent tests: skipped automatically when the atlas is not present
# ---------------------------------------------------------------------------


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
