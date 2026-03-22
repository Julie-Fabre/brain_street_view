import os
import json
import requests


def fetch_connectivity_summary(experiment_id, save_file_path):
    url = f'http://connectivity.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json?criteria=[section_data_set_id$eq{experiment_id}]&num_rows=all'

    try:
        resp = requests.get(url)
        resp.raise_for_status()
    except Exception as e:
        print(f'Warning: Failed to get data for ID {experiment_id}: {e}')
        return False

    try:
        data = resp.json()
    except Exception as e:
        print(f'Warning: Failed to parse JSON for ID {experiment_id}: {e}')
        return False

    if not data.get('success'):
        print(f'Warning: API returned failure for experiment ID {experiment_id}')
        return False

    # Filter for injection data only
    injection_info = [entry for entry in data.get('msg', []) if entry.get('is_injection')]
    if not injection_info:
        print(f'Warning: No injection data found for experiment ID {experiment_id}')
        return False

    # Add experiment ID to each entry
    for entry in injection_info:
        entry['experimentID'] = experiment_id

    os.makedirs(save_file_path, exist_ok=True)
    with open(os.path.join(save_file_path, 'injectionSummary_all.json'), 'w') as f:
        json.dump(injection_info, f)

    return True


def load_injection_summary(save_file_path):
    with open(os.path.join(save_file_path, 'injectionSummary_all.json'), 'r') as f:
        return json.load(f)
