import time
import requests


def find_connectivity_experiments(regions, mouse_line='', primary_injection=True):
    """Query the Allen API for connectivity experiments by injection region.

    Parameters
    ----------
    regions : list of str
        Brain region acronyms to search for (e.g. ``['VISp', 'VISl']``).
    mouse_line : str, optional
        Transgenic mouse line filter. ``''`` for all lines, ``'0'`` for
        wild-type only.
    primary_injection : bool, optional
        If True, only return experiments where the region is the primary
        injection site.

    Returns
    -------
    list of int
        Experiment IDs matching the query.
    """
    experiment_ids = []
    primary = 'true' if primary_injection else 'false'

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

        if result is None or not result.get('success'):
            msg = result.get('msg', 'Unknown error') if result else 'Request failed'
            print(f"Query failed!\n{msg}\nAt URL: {url}")
            continue

        for entry in result.get('msg', []):
            if 'id' in entry:
                experiment_ids.append(entry['id'])

        print(f"Found {len(result.get('msg', []))} experiments in {region}")

    return experiment_ids
