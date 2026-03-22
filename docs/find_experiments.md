# Find connectivity experiments

Query the Allen API for anterograde tracing experiments by injection region, transgenic mouse line, and injection criteria.

## Python

```python
import bsv

# Find all experiments with injections in primary visual cortex
experiment_ids = bsv.find_connectivity_experiments(['VISp'])

# Find experiments across multiple visual areas
experiment_ids = bsv.find_connectivity_experiments(
    ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor'])

# Filter by transgenic mouse line ('' = all lines, '0' = wild-type)
experiment_ids = bsv.find_connectivity_experiments(['VISp'], mouse_line='0')

# Only primary injection sites
experiment_ids = bsv.find_connectivity_experiments(['VISp'], primary_injection=True)
```

## MATLAB

```matlab
% Find all experiments with injections in primary visual cortex
experimentIDs = bsv.findConnectivityExperiments({'VISp'});

% Find experiments across multiple visual areas
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl', 'VISal', 'VISam'});

% Filter by mouse line
experimentIDs = bsv.findConnectivityExperiments({'VISp'}, '0', true);
```

## Fetch and cache data

Once you have experiment IDs, download the projection density maps. Data is cached locally so subsequent calls skip the download.

### Python

```python
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, '/path/to/cache', 'my_query',
    'injectionIntensity', False,
    allen_atlas_path='/path/to/allenCCF')
```

### MATLAB

```matlab
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, 'my_query', 'injectionIntensity', false, '', allenAtlasPath);
```

### Normalization options

| Method | Description |
|---|---|
| `'none'` | Raw projection density values |
| `'injectionIntensity'` | Normalize by injection volume/intensity |

### Parameters

| Parameter | Description |
|---|---|
| `experiment_ids` | List of Allen experiment IDs |
| `save_location` | Local directory for caching downloaded data |
| `file_name` | Name for cached file (empty string = don't cache) |
| `normalization_method` | `'none'` or `'injectionIntensity'` |
| `subtract_other_hemisphere` | Subtract contralateral signal |
| `allen_atlas_path` | Path to local Allen CCF atlas files |
