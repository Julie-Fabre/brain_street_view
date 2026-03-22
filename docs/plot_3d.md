# Plot 3D projections

Render injection sites and target regions as 3D isosurfaces overlaid on atlas anatomy.

```{image} ../images/docs/plot_3d_VISam_CP.png
:width: 60%
:align: center
```
*3D rendering of VISam injection sites (dots) and caudate putamen (transparent volume).*

**Python:**
```python
bsv.plot_connectivity_3d(
    injection_summary=injection_summary,
    allen_atlas_path='/path/to/allenCCF',
    region_to_plot='CP',
    color=[[0.543, 0, 0], [0, 0.746, 1]],
    plot_patch=True,
    animate=True)
```

**MATLAB:**
```matlab
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, 'CP', color, true);
```
