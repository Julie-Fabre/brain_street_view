# Find and fetch connectivity data

Query the Allen API for anterograde tracing experiments, then download and cache the projection density maps.

**Python:**
```python
import bsv

# Find experiments with injections in visual cortex
experiment_ids = bsv.find_connectivity_experiments(
    regions=['VISp', 'VISl', 'VISal', 'VISam'],  # list of injection region acronyms to search for
    mouse_line='',           # '' = all lines, '0' = wild-type only
    primary_injection=True)  # only primary injection sites

# Download and cache projection data
(experiment_imgs,       # ndarray (AP x DV x ML x groups): averaged projection fluorescence volumes
 injection_summary,     # dict of lists: injection metadata per experiment (coordinates, volumes, etc.)
 _,                     # ndarray or None: per-experiment volumes (only if load_all=True)
 _                      # dict: per-experiment region and group metadata
) = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location='/path/to/cache',      # local directory for caching downloaded data
    file_name='my_query',                # base name for the cached metadata CSV ('' to skip)
    normalization_method='injectionIntensity',  # 'none' or 'injectionIntensity' (divide by injection volume)
    subtract_other_hemisphere=False,            # subtract contralateral hemisphere signal
    allen_atlas_path='/path/to/allenCCF')       # atlas files auto-downloaded here on first use
```

**MATLAB:**
```matlab
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl', 'VISal', 'VISam'}, '', true);
%                                                 injection region acronyms    mouse line  primary only

[experimentImgs, ...         % ndarray (AP x DV x ML x groups): averaged projection fluorescence volumes
 injectionSummary, ...       % struct: injection metadata per experiment (coordinates, volumes, etc.)
 individualProjections, ...  % ndarray or []: per-experiment volumes (only if loadAll=true)
 experimentRegionInfo ...    % struct: per-experiment region and group metadata
] = bsv.fetchConnectivityData( ...
    experimentIDs, ...       % experiment IDs from findConnectivityExperiments
    saveLocation, ...        % local directory for caching downloaded data
    'my_query', ...          % base name for cached metadata CSV ('' to skip)
    'injectionIntensity', ...% normalization: 'none' or 'injectionIntensity'
    false, ...               % subtract contralateral hemisphere signal
    '', ...                  % grouping method: 'AP', 'ML', 'DV', or '' for none
    allenAtlasPath);         % atlas files auto-downloaded here on first use
```

Data is cached locally after the first download, so subsequent calls load from disk.
