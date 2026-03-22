# Compare regions

Plot projections from grouped source regions, or to multiple target regions side by side.

```{image} ../images/docs/plot_grouped_VIS_CP.png
:width: 100%
```
*Projections to caudate putamen from grouped visual areas: VISp (top), VISl+VISal (middle), remaining visual areas (bottom).*

**Python:**
```python
# Group source regions and fetch with grouping
grouped_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
region_groups = [1, 2, 2, 3, 3, 3, 3, 3]

experiment_ids = bsv.find_connectivity_experiments(grouped_regions)
experiment_imgs, _, _, region_info = bsv.fetch_connectivity_data(
    experiment_ids, '/path/to/cache', '',
    'injectionIntensity', False,
    allen_atlas_path='/path/to/allenCCF',
    input_regions=grouped_regions, region_groups=region_groups)

# Plot — one row per group
bsv.plot_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity',
    input_regions=grouped_regions, region_groups=region_groups,
    experiment_region_info=region_info)

# Compare multiple target regions side by side
bsv.plot_connectivity_multi_region(
    experiment_imgs, '/path/to/allenCCF', ['CP', 'ACB'],
    10, 15, 'coronal', True, 2, 'global', color, 'injectionIntensity')
```

**MATLAB:**
```matlab
[experimentImgs, ~, ~, regionInfo] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, '', 'injectionIntensity', false, '', allenAtlasPath, false, ...
    {'VISp', 'VISl', 'VISal', 'VISam'}, [1, 1, 2, 2]);

bsv.plotConnectivityMultiRegion(experimentImgs, allenAtlasPath, ...
    {'CP', 'ACB'}, 10, 15, 'coronal', true, 2, 'global', color, 'injectionIntensity')
```
