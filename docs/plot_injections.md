# Plot injection sites

Visualize where viral injections were placed, either one region at a time or all regions together.

```{image} ../images/docs/plot_injections_VISam.png
:width: 100%
```
*Injection sites in antero-medial visual cortex (VISam), shown across coronal slices.*

## Plot injections by region

Show injection sites for each source region separately.

### Python

```python
import bsv

# Plot each input region's injections one by one
bsv.plot_multi_region_injections(
    experiment_imgs, '/path/to/allenCCF',
    ['VISp', 'VISl', 'VISal'],
    10, 15, 'coronal', True, 2, 'global',
    color, 'injectionIntensity')

# Or manually loop through regions
for region in ['VISp', 'VISl', 'VISal']:
    bsv.plot_connectivity(
        experiment_imgs, '/path/to/allenCCF', region,
        10, 15, 'coronal', True, 2, 'global',
        color, 'injectionIntensity')
```

### MATLAB

```matlab
% Plot each input region's injections
for iRegion = 1:length(inputRegions)
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iRegion), ...
        10, 15, 'coronal', true, 2, 'global', color, 'injectionIntensity')
end
```

## Combined injection view

Show injections from multiple source regions overlaid.

### Python

```python
bsv.plot_injections_combined(
    experiment_imgs, '/path/to/allenCCF',
    input_regions, color)
```

### MATLAB

```matlab
bsv.plotInjectionsCombined(experimentImgs, allenAtlasPath, inputRegions, color)
```
