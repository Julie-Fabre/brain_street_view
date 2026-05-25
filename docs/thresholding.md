# Threshold projections

Identify significant projection signals using absolute, percentile, z-score, or relative thresholds.

```{image} ../images/docs/threshold_VISam_CP.png
:width: 100%
```
*VISam projections to CP thresholded at the 90th percentile. Red dots mark voxels above threshold.*

**Python:**
```python
# Percentile - keep top 10% of voxels
bsv.threshold_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',  # atlas files auto-downloaded here on first use
    input_region='CP',                     # target region acronym
    number_of_chunks=10,                   # number of evenly spaced slices to display
    number_of_pixels=15,                   # number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    plane='coronal',                       # 'coronal' or 'sagittal'
    region_only=True,                      # mask display to target region boundary
    smoothing=2,                           # Gaussian smoothing sigma in bins (0 for none)
    color_limits='global',                 # 'global', 'per_slice', or [min, max]
    color=color,                           # RGB colour(s) for region groups
    threshold=90,                          # threshold value (interpretation depends on threshold_method)
    threshold_method='percentile')         # 'percentile', 'absolute', 'zscore', or 'relative'

# threshold_method options:
#   'percentile' - keep voxels above the Nth percentile (threshold=90 → top 10%)
#   'absolute'   - raw density cutoff
#   'zscore'     - z-score cutoff
#   'relative'   - fraction of max (threshold between 0 and 1)
```

**MATLAB:**
```matlab
bsv.thresholdConnectivity( ...
    experimentImgs, ...    % projection fluorescence array from fetchConnectivityData
    allenAtlasPath, ...    % atlas files auto-downloaded here on first use
    'CP', ...              % target region acronym
    10, ...                % number of evenly spaced slices to display
    15, ...                % number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    'coronal', ...         % 'coronal' or 'sagittal'
    true, ...              % region_only: mask to target region boundary
    2, ...                 % smoothing sigma in bins
    'global', ...          % color_limits: 'global', 'per_slice', or [min, max]
    color, ...             % RGB colour(s) for region groups
    90, ...                % threshold value
    'percentile')          % threshold_method: 'percentile', 'absolute', 'zscore', or 'relative'
```
