# Brain Street View

Load, visualize, and analyze [Allen Mouse Brain Connectivity Atlas](http://connectivity.brain-map.org/) data ([Oh et al., Nature, 2014](https://doi.org/10.1038/nature13186)).

Available in both **Python** and **MATLAB**.

📦 **Source code:** [github.com/Julie-Fabre/brain_street_view](https://github.com/Julie-Fabre/brain_street_view)

```{image} ../images/bsv.png
:width: 200px
:align: center
```

## Installation

**Python:** (3.9–3.12). We recommend a fresh conda environment:

```bash
conda create -n bsv python=3.11
conda activate bsv
pip install brain-street-view
```

A `venv` works equally well, and plain `pip install brain-street-view` is fine if you
don't want a dedicated environment.

Or from source:

```bash
conda create -n bsv python=3.11
conda activate bsv
git clone https://github.com/Julie-Fabre/brain_street_view.git
cd brain_street_view
pip install -e .
```

**MATLAB:** Clone the [repository](https://github.com/Julie-Fabre/brain_street_view) and add it to [MATLAB's path](https://uk.mathworks.com/help/matlab/ref/pathtool.html). Requires MATLAB >= 2019a. Also needs [npy-matlab](https://github.com/kwikteam/npy-matlab), [brewermap](https://github.com/DrosteEffect/BrewerMap), and [prettify-matlab](https://github.com/Julie-Fabre/prettify_matlab).

**Both:** The Allen CCF atlas files (~2.4 GB) are downloaded automatically the first time a plotting function is called and saved to your `allen_atlas_path`. You can also download them manually from [figshare](https://figshare.com/articles/dataset/Modified_Allen_CCF_2017_for_cortex-lab_allenCCF/25365829).

```{toctree}
:maxdepth: 1
:caption: Contents

getting_started
find_and_fetch
plot_2d
plot_3d
plot_injections
upstream_projectome
compare_regions
thresholding
striatum_analysis
examples
api
contributing
```
