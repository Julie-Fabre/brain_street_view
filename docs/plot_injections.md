# Plot injection sites

Visualize where viral injections were placed, either one region at a time or combined.

```{image} ../images/docs/plot_injections_VISam.png
:width: 100%
```
*Injection sites in antero-medial visual cortex (VISam), shown across coronal slices.*

**Python:**
```python
# Plot each region's injections separately (one row per region)
bsv.plot_multi_region_injections(
    experiment_imgs=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',       # atlas files auto-downloaded here on first use
    input_regions=['VISp', 'VISl', 'VISal'],    # source regions to display
    number_of_slices=10,                         # number of evenly spaced slices to display
    number_of_pixels=15,                         # number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    plane='coronal',                             # 'coronal' or 'sagittal'
    region_only=True,                            # mask display to source region boundary
    smoothing=2,                                 # Gaussian smoothing sigma in bins (0 for none)
    color_limits='global',                       # 'global', 'per_slice', or [min, max]
    color=color,                                 # RGB colour(s) per region group
    normalization_method='injectionIntensity')

# Or show all regions overlaid in a single row
bsv.plot_injections_combined(
    experiment_imgs=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    input_regions=input_regions,
    number_of_slices=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=color,
    normalization_method='injectionIntensity')
```

**MATLAB:**
```matlab
% One region at a time
for iRegion = 1:length(inputRegions)
    bsv.plotConnectivity( ...
        experimentImgs, ...         % projection fluorescence array from fetchConnectivityData
        allenAtlasPath, ...         % atlas files auto-downloaded here on first use
        inputRegions(iRegion), ...  % source region acronym to visualize
        10, ...                     % number of evenly spaced slices to display
        15, ...                     % number of 2D histogram bins per axis per slice (bin size adapts to region extent)
        'coronal', ...              % 'coronal' or 'sagittal'
        true, ...                   % region_only: mask to source region boundary
        2, ...                      % smoothing sigma in bins
        'global', ...               % color_limits
        color, ...                  % RGB colour(s) per region group
        'injectionIntensity')
end

% Combined (all regions overlaid)
bsv.plotInjectionsCombined(experimentImgs, allenAtlasPath, inputRegions, color)
```
