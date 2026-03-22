# Find and fetch connectivity data

Query the Allen API for anterograde tracing experiments, then download and cache the projection density maps.

**Python:**
```python
import bsv

# Find experiments with injections in visual cortex
experiment_ids = bsv.find_connectivity_experiments(
    regions=['VISp', 'VISl', 'VISal', 'VISam'],
    mouse_line='',           # '' = all lines, '0' = wild-type only
    primary_injection=True)  # only primary injection sites

# Download and cache projection data
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location='/path/to/cache',
    file_name='my_query',
    normalization_method='injectionIntensity',
    subtract_other_hemisphere=False,
    allen_atlas_path='/path/to/allenCCF')
```

**MATLAB:**
```matlab
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl', 'VISal', 'VISam'}, '', true);

[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, 'my_query', 'injectionIntensity', false, '', allenAtlasPath);
```

Data is cached locally after the first download, so subsequent calls load from disk.
