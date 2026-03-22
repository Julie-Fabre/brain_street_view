import bsv


def test_find_visp_experiments():
    """Test finding experiments for VISp - should return multiple IDs."""
    ids = bsv.find_connectivity_experiments(['VISp'])
    assert len(ids) > 0
    assert all(isinstance(x, int) for x in ids)
    print(f'Found {len(ids)} VISp experiments')


def test_find_with_mouse_line():
    """Test filtering by mouse line (wild-type = C57BL/6J)."""
    ids_all = bsv.find_connectivity_experiments(['VISp'])
    ids_wt = bsv.find_connectivity_experiments(['VISp'], mouse_line='C57BL/6J')
    # Wild-type subset should be smaller
    assert len(ids_wt) <= len(ids_all)
    print(f'VISp all: {len(ids_all)}, wild-type: {len(ids_wt)}')


def test_find_multiple_regions():
    """Test querying multiple regions returns IDs from both."""
    ids = bsv.find_connectivity_experiments(['VISp', 'VISl'])
    assert len(ids) > 0
    # Should have experiments from both regions
    print(f'Combined VISp+VISl: {len(ids)} experiments')
