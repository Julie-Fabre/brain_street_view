import time
import requests

from .atlas_utils import load_projection_info


def fetch_upstream_regions(target_regions, min_projection_volume=0.01, mouse_line=''):
    """Query the Allen API for brain regions that project TO given target region(s).

    Uses the precomputed ProjectionStructureUnionize table to find all anterograde
    tracing experiments with signal in the target structure(s), then returns the
    injection site acronyms (= upstream source regions).

    Parameters
    ----------
    target_regions : str or list of str
        Target brain region acronym(s) (e.g. ``'CP'`` or ``['CP', 'STR']``).
    min_projection_volume : float, optional
        Minimum projection_volume (mm^3) threshold to count as a projection.
        Default is 0.01.
    mouse_line : str, optional
        Reserved for API consistency with :func:`find_connectivity_experiments`.
        Not used in the PSU query.

    Returns
    -------
    list of str
        Unique injection-site acronyms of experiments projecting to the target,
        sorted alphabetically.
    """
    if isinstance(target_regions, str):
        target_regions = [target_regions]

    projection_info = load_projection_info()
    all_upstream = set()

    for target in target_regions:
        # 1. Resolve acronym -> Allen structure ID
        id_url = (
            'http://api.brain-map.org/api/v2/data/query.json'
            f"?criteria=model::Structure,rma::criteria,[acronym$eq'{target}']&num_rows=1"
        )
        structure_id = None
        for attempt in range(3):
            try:
                resp = requests.get(id_url)
                data = resp.json()
                if data.get('success') and data.get('msg'):
                    structure_id = data['msg'][0]['id']
                    break
            except Exception:
                pass
            if attempt < 2:
                time.sleep(1)

        if structure_id is None:
            print(f"Warning: Could not resolve structure ID for '{target}'")
            continue

        print(f"Querying upstream projections to {target} (structure_id={structure_id})...")

        # 2. Query ProjectionStructureUnionize for non-injection records in this structure
        psu_url = (
            'http://connectivity.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json'
            f'?criteria=[structure_id$eq{structure_id}][is_injection$eqfalse]'
            f'[projection_volume$gt{min_projection_volume}]'
            '&num_rows=all'
        )
        result = None
        for attempt in range(3):
            try:
                resp = requests.get(psu_url)
                result = resp.json()
                if result.get('success'):
                    break
            except Exception:
                pass
            if attempt < 2:
                time.sleep(1)

        if result is None or not result.get('success'):
            print(f"Warning: PSU query failed for target '{target}'")
            continue

        records = result.get('msg', [])
        exp_ids = {r['section_data_set_id'] for r in records}
        print(f"  Found {len(exp_ids)} experiments projecting to {target}")

        # 3. Look up injection structure acronyms from local projection_info CSV
        matched = projection_info[projection_info['id'].isin(exp_ids)]
        acronyms = matched['structure_abbrev'].dropna().unique().tolist()
        all_upstream.update(acronyms)

    return sorted(all_upstream)
