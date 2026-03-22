# Quick start

## Python

```python
import bsv

# 1. Find experiments — query Allen API for injection experiments
experiment_ids = bsv.find_connectivity_experiments(['VISp', 'VISl'])

# 2. Fetch data — download and cache projection density maps
imgs, inj_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, '/path/to/cache', '',
    'injectionIntensity', False,
    allen_atlas_path='/path/to/allenCCF')

# 3. Plot 2D projections to caudate putamen
bsv.plot_connectivity(imgs, '/path/to/allenCCF', 'CP',
                       10, 15, 'coronal', True, 2, 'global', None,
                       'injectionIntensity')

# 4. 3D visualization
bsv.plot_connectivity_3d(inj_summary, '/path/to/allenCCF', 'CP',
                          plot_patch=True)
```

See `example.ipynb` or `example.py` for the full workflow including region grouping, thresholding, and CP subregion analysis. The first run downloads images from the Allen API and caches them locally; subsequent runs load from cache.

## MATLAB

Requires MATLAB >= 2019a. Clone the repo and add it (plus dependencies) to MATLAB's path.

```matlab
% 1. Find experiments
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl'});

% 2. Fetch data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, fileName, 'injectionIntensity', false, '', allenAtlasPath);

% 3. Plot
bsv.plotConnectivity(experimentImgs, allenAtlasPath, 'CP', 10, 15, ...
    'coronal', true, 2, 'global', [], 'injectionIntensity')
```

See `+bsv/example.m` and `gettingStarted.mlx` for the full MATLAB workflow.

## Pipeline overview

Brain Street View follows a four-step pipeline:

1. **Find experiments** — Query the Allen API for connectivity experiments filtered by injection region, transgenic mouse line, and injection quality.
2. **Fetch data** — Download projection density maps, cache them locally, and apply normalization (injection intensity, per-region, z-score, or robust scaling).
3. **Visualize** — Generate 2D coronal/sagittal slice views or 3D isosurface renderings, with region masking, thresholding, and customizable colormaps.
4. **Analyze** — Break down projections by anatomical subdivisions (e.g., CP subregions), with per-slice and summary statistics and optional CSV export.
