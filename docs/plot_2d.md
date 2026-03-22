# Plot 2D projections

Visualize projection density as 2D coronal or sagittal slices, masked to a target region.

```{image} ../images/docs/plot_2d_VISam_CP.png
:width: 100%
```
*Projections from antero-medial visual cortex (VISam) to caudate putamen (CP), shown across 10 evenly spaced coronal slices.*

**Python:**
```python
proj_array, proj_coords = bsv.plot_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10,           # number of slices
    15,           # pixels per slice
    'coronal',    # 'coronal' or 'sagittal'
    True,         # mask to target region only
    2,            # smoothing in pixels (0 = none)
    'global',     # color scale: 'global', 'per_slice', or [min, max]
    None,         # color (None = default)
    'injectionIntensity')
```

**MATLAB:**
```matlab
[projArray, projCoords] = bsv.plotConnectivity(experimentImgs, allenAtlasPath, ...
    'CP', 10, 15, 'coronal', true, 2, 'global', [], 'injectionIntensity');
```

Use `custom_slices` to plot specific slice indices instead of evenly spaced ones.
