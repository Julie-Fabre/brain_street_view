import time
import requests

from .fetch_upstream_regions import fetch_upstream_regions
from .atlas_utils import load_projection_info


def find_connectivity_experiments(regions, mouse_line='', primary_injection=True, target_regions=None):
    """Query the Allen API for connectivity experiments by injection region.

    Parameters
    ----------
    regions : list of str
        Brain region acronyms to search for (e.g. ``['VISp', 'VISl']``).
        Ignored when *target_regions* is provided and *regions* is empty.
    mouse_line : str, optional
        Transgenic mouse line filter. ``''`` for all lines, ``'0'`` for
        wild-type only.
    primary_injection : bool, optional
        If True, only return experiments where the region is the primary
        injection site.
    target_regions : str or list of str or None, optional
        If provided, call :func:`fetch_upstream_regions` to discover which
        source regions project to these target(s), then query experiments
        injected in those sources. If *regions* is also non-empty it is used
        as an additional filter (intersection).

    Returns
    -------
    list of int
        Experiment IDs matching the query.
    """
    if target_regions is not None:
        upstream = fetch_upstream_regions(target_regions, mouse_line=mouse_line)
        if regions:
            # Match upstream leaves against caller-supplied regions, supporting
            # parent acronyms: 'ACA' should match 'ACAd', 'ACAv', etc.
            upstream = [
                r for r in upstream
                if any(r == s or r.startswith(s) for s in regions)
            ]
        regions = upstream
        if not regions:
            print("Warning: No upstream regions found for the specified target(s).")
            return []

    experiment_ids = []
    primary = 'true' if primary_injection else 'false'
    projection_info = load_projection_info()

    # Expand parent regions to their children.
    # A region is considered a parent if it never appears as a primary injection
    # site in the database but prefix-matches child regions that do (same
    # convention as find_structure_indices in atlas_utils).
    all_abbrevs = set(projection_info['structure_abbrev'].dropna())
    expanded = []
    for region in regions:
        if region in all_abbrevs:
            expanded.append(region)
        else:
            children = sorted(
                a for a in all_abbrevs
                if a != region and a.startswith(region)
            )
            if children:
                print(f"'{region}' expanded to children: {children}")
                expanded.extend(children)
            else:
                expanded.append(region)  # pass through; API will handle it
    # Deduplicate while preserving order
    seen = set()
    regions = [r for r in expanded if not (r in seen or seen.add(r))]

    for region in regions:
        url = 'http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure'
        url += f'[injection_structures$eq{region}]'
        if mouse_line:
            url += f'[transgenic_lines$eq{mouse_line}]'
        url += f'[primary_structure_only$eq{primary}]'

        result = None
        for attempt in range(3):
            try:
                resp = requests.get(url)
                result = resp.json()
                if result.get('success'):
                    break
            except Exception:
                pass
            if attempt < 2:
                time.sleep(1)

        if result is not None and result.get('success'):
            for entry in result.get('msg', []):
                if 'id' in entry:
                    experiment_ids.append(entry['id'])
            print(f"Found {len(result.get('msg', []))} experiments in {region}")
        else:
            # Allen service fails for some valid regions; fall back to local CSV
            rows = projection_info[projection_info['structure_abbrev'] == region]
            if mouse_line:
                rows = rows[rows['transgenic_line'] == mouse_line]
            fallback_ids = rows['id'].tolist()
            experiment_ids.extend(fallback_ids)
            print(f"Found {len(fallback_ids)} experiments in {region} (via local CSV fallback)")

    return experiment_ids
