# Plot injection sites

Visualize where viral injections were placed, either one region at a time or combined.

```{image} ../images/docs/plot_injections_VISam.png
:width: 100%
```
*Injection sites in antero-medial visual cortex (VISam), shown across coronal slices.*

**Python:**
```python
# Plot each region's injections separately
bsv.plot_multi_region_injections(
    experiment_imgs=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    input_regions=['VISp', 'VISl', 'VISal'],
    number_of_slices=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=color,
    normalization_method='injectionIntensity')

# Or show all regions overlaid
bsv.plot_injections_combined(
    experiment_imgs=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    input_regions=input_regions,
    color=color)
```

**MATLAB:**
```matlab
% One region at a time
for iRegion = 1:length(inputRegions)
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iRegion), ...
        10, 15, 'coronal', true, 2, 'global', color, 'injectionIntensity')
end

% Combined
bsv.plotInjectionsCombined(experimentImgs, allenAtlasPath, inputRegions, color)
```
