"""Network-free tests for fetch_upstream_regions and find_connectivity_experiments.

The Allen API HTTP calls are mocked, so these run deterministically in CI without
network access or the atlas. They exercise the parsing/joining/fallback logic that
sits around the API calls.

Note: each bsv submodule is named after the function it exports (e.g. the function
``bsv.fetch_upstream_regions`` shadows the submodule of the same name in the package
namespace). We therefore fetch the real module objects from ``sys.modules`` and patch
them with ``patch.object`` rather than a dotted string target.
"""
import sys
from unittest.mock import patch, MagicMock

import bsv
from bsv.fetch_upstream_regions import fetch_upstream_regions

# Real module objects (not the shadowing functions of the same name).
_FUR = sys.modules['bsv.fetch_upstream_regions']
_FCE = sys.modules['bsv.find_connectivity_experiments']


def _json_response(payload):
    """A fake requests.Response whose .json() returns *payload*."""
    resp = MagicMock()
    resp.json.return_value = payload
    return resp


# Two experiment IDs that are VISp injections in the bundled projection_info CSV.
_VISP_EXPERIMENT_IDS = [180296424, 114008926]


def test_fetch_upstream_regions_joins_psu_to_local_csv():
    """PSU records are joined to the local CSV to recover injection acronyms."""
    structure_resp = _json_response({'success': True, 'msg': [{'id': 672}]})  # CP id
    psu_resp = _json_response({
        'success': True,
        'msg': [{'section_data_set_id': i} for i in _VISP_EXPERIMENT_IDS],
    })

    # First GET resolves the structure id, second is the PSU query.
    with patch.object(_FUR, 'requests') as mock_requests:
        mock_requests.get.side_effect = [structure_resp, psu_resp]
        upstream = fetch_upstream_regions('CP')

    # Both IDs map to VISp injections -> deduplicated, sorted single entry.
    assert upstream == ['VISp']


def test_fetch_upstream_regions_accepts_string_or_list():
    """A bare string target is coerced to a one-element list (no crash)."""
    structure_resp = _json_response({'success': True, 'msg': [{'id': 672}]})
    psu_resp = _json_response({'success': True, 'msg': []})  # no projections found
    with patch.object(_FUR, 'requests') as mock_requests:
        mock_requests.get.side_effect = [structure_resp, psu_resp]
        assert fetch_upstream_regions('CP') == []


def test_fetch_upstream_regions_unresolved_structure_returns_empty():
    """When the structure id cannot be resolved the function degrades gracefully."""
    empty_resp = _json_response({'success': True, 'msg': []})
    with patch.object(_FUR, 'requests') as mock_requests, \
            patch.object(_FUR, 'time'):  # skip retry backoff
        mock_requests.get.return_value = empty_resp
        assert fetch_upstream_regions('NOTAREGION') == []


def test_find_connectivity_experiments_empty_upstream_returns_empty():
    """target_regions with no upstream sources yields an empty result, not an error."""
    with patch.object(_FCE, 'fetch_upstream_regions', return_value=[]):
        assert bsv.find_connectivity_experiments([], target_regions='CP') == []


def test_find_connectivity_experiments_falls_back_to_local_csv():
    """If the Allen service query fails, IDs are recovered from the bundled CSV."""
    failed = _json_response({'success': False})
    with patch.object(_FCE, 'requests') as mock_requests, \
            patch.object(_FCE, 'time'):
        mock_requests.get.return_value = failed
        ids = bsv.find_connectivity_experiments(['VISp'])

    assert len(ids) > 0
    assert all(isinstance(x, int) for x in ids)
