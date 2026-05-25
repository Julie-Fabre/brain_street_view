
# Brain Street View <img src="./images/bsv.png" width="20%" title="bsv" alt="bsv" align="left" vspace = "20">

[![PyPI version](https://img.shields.io/pypi/v/brain-street-view.svg)](https://pypi.org/project/brain-street-view/)
[![Tests](https://github.com/Julie-Fabre/brain_street_view/actions/workflows/tests.yml/badge.svg)](https://github.com/Julie-Fabre/brain_street_view/actions/workflows/tests.yml)
[![Documentation](https://readthedocs.org/projects/brain-street-view/badge/?version=latest)](https://brain-street-view.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![MATLAB 2019a+](https://img.shields.io/badge/MATLAB-2019a%2B-orange.svg)](https://mathworks.com/products/matlab.html)

Load and plot Allen Connectivity Data ([Oh et al., Nature, 2014](https://doi.org/10.1038/nature13186))

Available in both **MATLAB** and **Python**. **[Documentation](https://brain-street-view.readthedocs.io)**

### 🏁 Quick start

**Python**:
```bash
pip install brain-street-view
```
Then open `example.ipynb` or run `example.py`. The first run downloads images from the Allen API and caches them locally. Subsequent runs load from cache.

**MATLAB**:
See the script `gettingStarted.mlx`. Requires MATLAB>=2019a.

### ⚒️ Installation

#### Python

```bash
pip install brain-street-view
```

Or from source:
```bash
git clone https://github.com/Julie-Fabre/brain_street_view.git
cd brain_street_view
pip install -e .
```

You also need the Allen CCF atlas files (not included):
- [allenCCF](https://github.com/cortex-lab/allenCCF) — annotation volumes and structure trees

#### MATLAB

- [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the [repository](https://github.com/Julie-Fabre/brain_street_view) and the dependencies below.
- Add BrainStreetView's and the dependencies' folders to [MATLAB's path](https://uk.mathworks.com/help/matlab/ref/pathtool.html).

MATLAB dependencies:
- [allenCCF](https://github.com/cortex-lab/allenCCF), to get Allen Atlas files
- [npy-matlab](https://github.com/kwikteam/npy-matlab), to read in .npy files
- [brewermap](https://github.com/DrosteEffect/BrewerMap), to generate colormaps
- [prettify-matlab](https://github.com/Julie-Fabre/prettify_matlab), to make plots pretty.

### 📖 Usage

#### Python
```python
import bsv

# 1. Find experiments by injected region (source search)
# Checkout region acronyms here: https://connectivity.brain-map.org/
experiment_ids = bsv.find_connectivity_experiments(
    regions=['VISp', 'VISl'])       # list of injection region acronyms to search for

# 2. Fetch fluorescence data
(imgs,                       # ndarray (AP x DV x ML x groups): averaged projection fluorescence volumes
 inj_summary,               # dict of lists: injection metadata per experiment (coordinates, volumes, etc.)
 individual_projections,    # ndarray or None: per-experiment volumes (only populated if load_all=True)
 experiment_region_info     # dict: per-experiment region and group metadata
) = bsv.fetch_connectivity_data(
    experiment_ids=experiment_ids,
    save_location='/path/to/cache', # local directory for caching downloaded data
    file_name='',                   # base name for the cached metadata CSV ('' to skip)
    normalization_method='injectionIntensity',  # 'none' or 'injectionIntensity' (divide by injection volume)
    subtract_other_hemisphere=False,            # subtract contralateral hemisphere signal
    allen_atlas_path='/path/to/allenCCF')       # path to the Allen CCF atlas directory

# 3. Plot projections to striatum
bsv.plot_connectivity(
    experiment_data=imgs,
    allen_atlas_path='/path/to/allenCCF',  # path to the Allen CCF atlas directory
    output_region='CP',                    # target region acronym to visualize
    number_of_chunks=10,                   # number of evenly spaced slices to display
    number_of_pixels=15,                   # pixel resolution per slice panel
    plane='coronal',                       # 'coronal' or 'sagittal'
    region_only=True,                      # mask display to target region boundary
    smoothing=2,                           # Gaussian smoothing sigma in pixels (0 for none)
    color_limits='global',                 # 'global', 'per_slice', or [min, max]
    color=None,                            # RGB colour(s) for region groups, or None for default
    normalization_info='injectionIntensity')

# 4. 3D visualization
bsv.plot_connectivity_3d(
    injection_summary=inj_summary,
    allen_atlas_path='/path/to/allenCCF',  # path to the Allen CCF atlas directory
    region_to_plot='CP',                   # target region acronym
    plot_patch=True)                       # render region as solid isosurface (False for grid)
```

See `example.ipynb` for the full workflow including region grouping, thresholding, and CP subregion analysis.

#### MATLAB
```matlab
% 1. Find experiments
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl'}); % list of injection region acronyms

% 2. Fetch data
[experimentImgs, ...         % ndarray (AP x DV x ML x groups): averaged projection fluorescence volumes
 injectionSummary, ...      % struct: injection metadata per experiment (coordinates, volumes, etc.)
 individualProjections, ... % ndarray or []: per-experiment volumes (only populated if loadAll=true)
 experimentRegionInfo ...   % struct: per-experiment region and group metadata
] = bsv.fetchConnectivityData( ...
    experimentIDs, ...      % experiment IDs from findConnectivityExperiments
    saveLocation, ...       % local directory for caching downloaded data
    fileName, ...           % base name for cached metadata CSV ('' to skip)
    'injectionIntensity', ...  % normalization: 'none' or 'injectionIntensity'
    false, ...              % subtract contralateral hemisphere signal
    '', ...                 % grouping method: 'AP', 'ML', 'DV', or '' for none
    allenAtlasPath);        % path to the Allen CCF atlas directory

% 3. Plot
bsv.plotConnectivity( ...
    experimentImgs, ...     % projection density array from fetchConnectivityData
    allenAtlasPath, ...     % path to the Allen CCF atlas directory
    'CP', ...               % target region acronym to visualize
    10, ...                 % number of evenly spaced slices to display
    15, ...                 % pixel resolution per slice panel
    'coronal', ...          % plane: 'coronal' or 'sagittal'
    true, ...               % region_only: mask display to target region boundary
    2, ...                  % smoothing: Gaussian sigma in pixels (0 for none)
    'global', ...           % color_limits: 'global', 'per_slice', or [min, max]
    [], ...                 % color: RGB colour(s) for region groups, or [] for default
    'injectionIntensity')   % normalization method used during fetch
```

See `+bsv/example.m` for the full MATLAB workflow.

### 🖼️ Gallery
 - projections from visual cortices to striatum: injection sites and striatum plotted in 3D
<img src="./images/visualInjections.gif" width="50%" title="ex_vis_cp" vspace = "20">

- example injections in primary visual cortex (VISp)
<img src="./images/ex_VISp.png" width="100%" title="ex_visp" vspace = "20">

- example injections in antero-lateral visual cortex (VISal)
<img src="./images/ex_VISal.png" width="100%" title="ex_visal" vspace = "20">

- projections from visual cortices (VIS) to striatum (CP)
<img src="./images/ex_VIS_to_CP.png" width="100%" title="ex_vis_cp" vspace = "20">


### 📬 Contact me
If you run into any issues or if you have any suggestions, please raise a github issue, create a pull request or email me: [juliemfabre[at]gmail[dot]com](mailto:juliemfabre@gmail.com).
