# Plot injection sites

Visualize where viral injections were placed, either one region at a time or all regions together.

```{image} ../images/ex_VISp.png
:width: 100%
```
*Example injections in primary visual cortex (VISp).*

```{image} ../images/ex_VISal.png
:width: 100%
```
*Example injections in antero-lateral visual cortex (VISal).*

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
