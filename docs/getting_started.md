# Getting Started

This tutorial walks through a complete Brain Street View workflow: finding
experiments, downloading data, and producing publication-ready figures of
projection patterns to the striatum.

## Prerequisites

Install the package:

```bash
pip install brain-street-view
```

You also need the [Allen CCF atlas files](https://github.com/cortex-lab/allenCCF)
downloaded locally. Set the path once and reuse it:

```python
allen_atlas_path = '/path/to/allenCCF'
save_location = '/path/to/cache'   # where downloaded data will be stored
```

## Step 1: Find experiments

Query the Allen API for anterograde tracing experiments injected in one or
more brain regions. Here we search for injections in primary visual cortex
(VISp) and lateral visual cortex (VISl):

```python
import bsv

experiment_ids = bsv.find_connectivity_experiments(
    regions=['VISp', 'VISl'],
    mouse_line='',            # '' = all transgenic lines
    primary_injection=True)   # only primary injection sites

print(f'Found {len(experiment_ids)} experiments')
```

To restrict to wild-type mice only, pass ``mouse_line='0'``.

## Step 2: Download and cache projection data

Fetch the projection density volumes for every experiment. Data is
downloaded once from the Allen API and cached locally for future runs:

```python
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location=save_location,
    file_name='visual_cortex_query',
    normalization_method='injectionIntensity',
    subtract_other_hemisphere=False,
    allen_atlas_path=allen_atlas_path)
```

The returned ``experiment_imgs`` array contains the averaged projection
density across experiments (shape: AP x DV x ML). ``injection_summary``
holds per-experiment injection metadata used by the 3D viewer.

## Step 3: Visualize projections in 2D

Plot the mean projection density to the caudate putamen (CP) across 10
evenly spaced coronal slices:

```python
proj_array, proj_coords = bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path=allen_atlas_path,
    output_region='CP',
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    normalization_info='injectionIntensity')
```

Each panel shows the mean projection density within the CP boundary, with
darker values indicating stronger projections. Change ``plane='sagittal'``
for a sagittal view.

## Step 4: 3D rendering

Overlay injection sites on a transparent isosurface of the target region:

```python
bsv.plot_connectivity_3d(
    injection_summary=injection_summary,
    allen_atlas_path=allen_atlas_path,
    region_to_plot='CP',
    plot_patch=True,
    animate=True)
```

In a Jupyter notebook this returns a rotating HTML animation. Set
``animate=False`` to get a static matplotlib figure instead.

## Step 5: Threshold significant projections

Identify voxels with strong projection signal using a percentile threshold:

```python
bsv.threshold_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path=allen_atlas_path,
    input_region='CP',
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    threshold=90,
    threshold_method='percentile')
```

Other methods: ``'absolute'`` (raw density cutoff), ``'zscore'``, or
``'relative'`` (fraction of maximum).

## Step 6: Compare grouped regions

Group source regions and visualise one row per group:

```python
grouped_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm']
region_groups = [1, 2, 2, 3, 3, 3]

experiment_ids = bsv.find_connectivity_experiments(regions=grouped_regions)
experiment_imgs, _, _, region_info = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location=save_location,
    file_name='',
    normalization_method='injectionIntensity',
    subtract_other_hemisphere=False,
    allen_atlas_path=allen_atlas_path,
    input_regions=grouped_regions,
    region_groups=region_groups)

bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path=allen_atlas_path,
    output_region='CP',
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    normalization_info='injectionIntensity',
    input_regions=grouped_regions,
    region_groups=region_groups,
    experiment_region_info=region_info)
```

## Next steps

- See {doc}`striatum_analysis` for CP/NAc subregion breakdown
- See {doc}`examples` for published research use-cases
- See {doc}`api` for the full API reference
