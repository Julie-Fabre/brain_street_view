# Plot 2D projections

Visualize projection density as 2D coronal or sagittal slices, masked to a target region.

```{image} ../images/docs/plot_2d_VISam_CP.png
:width: 100%
```
*Projections from antero-medial visual cortex (VISam) to caudate putamen (CP), shown across 10 evenly spaced coronal slices.*

**Python:**
```python
(proj_array,    # ndarray (n_slices x n_bins x n_bins x groups): binned projection density per slice
 proj_coords    # list of bin edge arrays: spatial coordinates of each slice panel
) = bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',  # atlas files auto-downloaded here on first use
    output_region='CP',                    # target region acronym to visualize
    number_of_chunks=10,                   # number of evenly spaced slices to display
    number_of_pixels=15,                   # number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    plane='coronal',                       # 'coronal' or 'sagittal'
    region_only=True,                      # mask display to target region boundary
    smoothing=2,                           # Gaussian smoothing sigma in bins (0 for none)
    color_limits='global',                 # 'global', 'per_slice', or [min, max]
    color=None,                            # RGB colour(s) for region groups, or None for default
    normalization_info='injectionIntensity')
```

**MATLAB:**
```matlab
[projArray, ...   % ndarray (n_slices x n_bins x n_bins x groups): binned projection density per slice
 projCoords ...   % cell array of bin edge arrays: spatial coordinates of each slice panel
] = bsv.plotConnectivity( ...
    experimentImgs, ...    % projection fluorescence array from fetchConnectivityData
    allenAtlasPath, ...    % atlas files auto-downloaded here on first use
    'CP', ...              % target region acronym to visualize
    10, ...                % number of evenly spaced slices to display
    15, ...                % number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    'coronal', ...         % 'coronal' or 'sagittal'
    true, ...              % region_only: mask display to target region boundary
    2, ...                 % smoothing: Gaussian sigma in bins (0 for none)
    'global', ...          % color_limits: 'global', 'per_slice', or [min, max]
    [], ...                % color: RGB colour(s) for region groups, or [] for default
    'injectionIntensity'); % normalization method used during fetch
```

Use `custom_slices` to plot specific slice indices instead of evenly spaced ones.
