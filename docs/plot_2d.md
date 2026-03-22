# Plot 2D projections

Visualize projection density as 2D coronal or sagittal slices, masked to a target region.

```{image} ../images/docs/plot_2d_VISam_CP.png
:width: 100%
```
*Projections from antero-medial visual cortex (VISam) to caudate putamen (CP), shown across 10 evenly spaced coronal slices.*

**Python:**
```python
proj_array, proj_coords = bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    output_region='CP',
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    normalization_info='injectionIntensity')
```

**MATLAB:**
```matlab
[projArray, projCoords] = bsv.plotConnectivity(experimentImgs, allenAtlasPath, ...
    'CP', 10, 15, 'coronal', true, 2, 'global', [], 'injectionIntensity');
```

Use `custom_slices` to plot specific slice indices instead of evenly spaced ones.
