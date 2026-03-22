# Compare regions

Plot projections to multiple target regions side by side.

## Python

```python
import bsv

# Compare projections to CP and NAc simultaneously
bsv.plot_connectivity_multi_region(
    experiment_imgs, '/path/to/allenCCF',
    ['CP', 'ACB'],            # multiple output regions
    10, 15, 'coronal', True, 2, 'global',
    color, 'injectionIntensity')
```

## MATLAB

```matlab
bsv.plotConnectivityMultiRegion(experimentImgs, allenAtlasPath, ...
    {'CP', 'ACB'}, 10, 15, 'coronal', true, 2, 'global', color, 'injectionIntensity')
```

## Region grouping

Group source regions and normalize within groups to compare projection patterns across inputs.

### Python

```python
# Group visual areas and normalize per group
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, '/path/to/cache', 'grouped_query',
    'injectionIntensity', False,
    grouping_method='perGroup',
    allen_atlas_path='/path/to/allenCCF',
    input_regions=['VISp', 'VISl', 'VISal', 'VISam'],
    region_groups=[1, 1, 2, 2])  # group VISp+VISl and VISal+VISam

# Then plot with group-aware normalization
bsv.plot_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    'injectionIntensity',
    input_regions=['VISp', 'VISl', 'VISal', 'VISam'],
    region_groups=[1, 1, 2, 2],
    normalize_by_group=True)
```

### MATLAB

```matlab
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, 'grouped_query', 'injectionIntensity', false, 'perGroup', ...
    allenAtlasPath, false, {'VISp', 'VISl', 'VISal', 'VISam'}, [1, 1, 2, 2]);
```
