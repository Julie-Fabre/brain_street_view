# Plot 3D projections

Render injection sites and projection targets as 3D isosurfaces overlaid on atlas anatomy.

```{image} ../images/docs/plot_3d_VISam_CP.png
:width: 60%
:align: center
```
*3D rendering of VISam injection sites (dots) and caudate putamen (transparent volume).*

## Python

```python
import bsv

# 3D isosurface rendering of projections to CP
bsv.plot_connectivity_3d(
    injection_summary, '/path/to/allenCCF', 'CP',
    color=[[0.543, 0, 0], [0, 0.746, 1]],
    plot_patch=True)

# Grid-based rendering (alternative to patch)
bsv.plot_connectivity_3d(
    injection_summary, '/path/to/allenCCF', 'CP',
    plot_patch=False)
```

## MATLAB

```matlab
% Patch-based 3D rendering
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, 'CP', ...
    color, true)

% Grid-based 3D rendering
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, 'CP', ...
    color, false)
```

## Options

| Parameter | Description |
|---|---|
| `plot_patch` | `true`: solid isosurface volume; `false`: grid point rendering |
| `color` | RGB color(s) per experiment group |
| `animate` | (Python only) Generate a rotating animation |
| `atlas_resolution` | Atlas resolution in microns (default: 10) |
