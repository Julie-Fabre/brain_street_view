# Plot 2D projections

Visualize projection density as 2D coronal or sagittal slices, masked to a target region of interest.

```{image} ../images/docs/plot_2d_VISam_CP.png
:width: 100%
```
*Projections from antero-medial visual cortex (VISam) to caudate putamen (CP), shown across 10 evenly spaced coronal slices.*

## Python

```python
import bsv

# Plot projections to caudate putamen, 10 evenly spaced coronal slices
bsv.plot_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10,           # number of slices
    15,           # pixels per slice
    'coronal',    # plane: 'coronal' or 'sagittal'
    True,         # region_only: mask to target region
    2,            # smoothing (pixels)
    'global',     # color_limits: 'global', 'per_slice', or [min, max]
    None,         # color (None = default)
    'injectionIntensity')
```

## MATLAB

```matlab
bsv.plotConnectivity(experimentImgs, allenAtlasPath, 'CP', ...
    10, 15, 'coronal', true, 2, 'global', [], 'injectionIntensity')
```

## Options

| Parameter | Values | Description |
|---|---|---|
| `plane` | `'coronal'`, `'sagittal'` | Slice orientation |
| `region_only` | `true` / `false` | Mask output to the target region |
| `color_limits` | `'global'`, `'per_slice'`, `[min, max]` | How to scale the colormap |
| `smoothing` | integer | Gaussian smoothing in pixels (0 = none) |
| `normalization_method` | `'none'`, `'injectionIntensity'` | Label for the normalization used during data fetch |
| `custom_slices` | list of slice indices | Plot specific slices instead of evenly spaced |
