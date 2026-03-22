
# Brain Street View <img src="./images/bsv.png" width="20%" title="bsv" alt="bsv" align="left" vspace = "20">
Load and plot Allen Connectivity Data ([Oh et al., Nature, 2014](https://doi.org/10.1038/nature13186))

Available in both **MATLAB** and **Python**.

### 🏁 Quick start

**Python** (recommended):
```bash
pip install -e .
```
Then open `example.ipynb` or run `example.py`. The first run downloads images from the Allen API and caches them locally. Subsequent runs load from cache.

**MATLAB**:
See the script `gettingStarted.mlx`. Requires MATLAB>=2019a.

### ⚒️ Installation

#### Python

```bash
# clone and install
git clone https://github.com/Julie-Fabre/brain_street_view.git
cd brain_street_view
pip install -e .
```

Or with conda:
```bash
conda create -n bsv python=3.12 -y
conda activate bsv
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

# 1. Find experiments
experiment_ids = bsv.find_connectivity_experiments(['VISp', 'VISl'])

# 2. Fetch data
imgs, inj_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, '/path/to/cache', '',
    'injectionIntensity', False,
    allen_atlas_path='/path/to/allenCCF')

# 3. Plot projections to striatum
bsv.plot_connectivity(imgs, '/path/to/allenCCF', 'CP',
                       10, 15, 'coronal', True, 2, 'global', None,
                       'injectionIntensity')

# 4. 3D visualization
bsv.plot_connectivity_3d(inj_summary, '/path/to/allenCCF', 'CP',
                          plot_patch=True)
```

See `example.ipynb` for the full workflow including region grouping, thresholding, and CP subregion analysis.

#### MATLAB
```matlab
% 1. Find experiments
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl'});

% 2. Fetch data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, fileName, 'injectionIntensity', false, '', allenAtlasPath);

% 3. Plot
bsv.plotConnectivity(experimentImgs, allenAtlasPath, 'CP', 10, 15, ...
    'coronal', true, 2, 'global', [], 'injectionIntensity')
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
If you run into any issues or if you have any suggestions, please raise a github issue, create a pull request or email me: [juliemfabre[at]gmail[dot]com](mailto:julie.mfabre@gmail.com).
