# Compare regions

Plot projections from grouped source regions, or to multiple target regions side by side.

```{image} ../images/docs/plot_grouped_VIS_CP.png
:width: 100%
```
*Projections to caudate putamen from grouped visual areas: VISp (top), VISl+VISal (middle), remaining visual areas (bottom).*

**Python:**
```python
# Group source regions (region_groups assigns each region to a display row)
grouped_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
region_groups   = [1,      2,      2,       3,       3,       3,       3,       3      ]
#                  row 1   row 2   row 2    row 3    row 3    row 3    row 3    row 3

experiment_ids = bsv.find_connectivity_experiments(regions=grouped_regions)
(experiment_imgs,   # ndarray (AP x DV x ML x n_groups): one volume per display group
 _,
 _,
 region_info        # dict: per-experiment region and group metadata (pass to plot_connectivity)
) = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location='/path/to/cache',
    file_name='',
    normalization_method='injectionIntensity',
    subtract_other_hemisphere=False,
    allen_atlas_path='/path/to/allenCCF',  # atlas files auto-downloaded here on first use
    input_regions=grouped_regions,
    region_groups=region_groups)

# Plot - one row per group
bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    output_region='CP',
    number_of_chunks=10,                   # number of evenly spaced slices to display
    number_of_pixels=15,                   # number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    normalization_info='injectionIntensity',
    input_regions=grouped_regions,
    region_groups=region_groups,
    experiment_region_info=region_info)    # pass through for correct group averaging

# Compare multiple target regions side by side
bsv.plot_connectivity_multi_region(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    output_regions=['CP', 'ACB'],          # list of target region acronyms
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=color,
    normalization_info='injectionIntensity')
```

**MATLAB:**
```matlab
[experimentImgs, ~, ~, regionInfo] = bsv.fetchConnectivityData( ...
    experimentIDs, saveLocation, '', 'injectionIntensity', false, '', allenAtlasPath, false, ...
    {'VISp', 'VISl', 'VISal', 'VISam'}, ...  % input regions
    [1, 1, 2, 2]);                            % group assignments (one row per unique value)

bsv.plotConnectivityMultiRegion( ...
    experimentImgs, allenAtlasPath, ...
    {'CP', 'ACB'}, ...  % list of target region acronyms
    10, ...             % number of evenly spaced slices
    15, ...             % number of 2D histogram bins per axis per slice (bin size adapts to region extent)
    'coronal', true, 2, 'global', color, 'injectionIntensity')
```
