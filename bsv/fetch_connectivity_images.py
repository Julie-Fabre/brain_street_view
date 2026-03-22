import os
import zipfile
import tempfile
import requests

PROJECTION_GRID_SIZE = (132, 80, 114)


def fetch_connectivity_images(experiment_id, save_file_path):
    """Download density.raw for an experiment. Returns True on success."""
    os.makedirs(save_file_path, exist_ok=True)

    raw_path = os.path.join(save_file_path, 'density.raw')
    if os.path.exists(raw_path):
        return True

    url = f'http://api.brain-map.org/grid_data/download/{experiment_id}?include=density'

    try:
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
    except Exception as e:
        print(f'Warning: Failed to download data for experiment {experiment_id}: {e}')
        return False

    # Save zip and extract
    zip_path = os.path.join(save_file_path, 'temp.zip')
    try:
        with open(zip_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(save_file_path)
    except Exception as e:
        print(f'Warning: Failed to extract data for experiment {experiment_id}: {e}')
        return False
    finally:
        if os.path.exists(zip_path):
            os.remove(zip_path)

    if not os.path.exists(raw_path):
        print(f'Warning: density.raw not found after extraction for experiment {experiment_id}')
        return False

    return True
