# Brain Street View

Load, visualize, and analyze [Allen Mouse Brain Connectivity Atlas](http://connectivity.brain-map.org/) data ([Oh et al., Nature, 2014](https://doi.org/10.1038/nature13186)).

Available in both **Python** and **MATLAB**.

```{image} ../images/bsv.png
:width: 200px
:align: center
```

## Installation

**Python:**

```bash
pip install brain-street-view
```

Or from source:

```bash
git clone https://github.com/Julie-Fabre/brain_street_view.git
cd brain_street_view
pip install -e .
```

**MATLAB:** Clone the [repository](https://github.com/Julie-Fabre/brain_street_view) and add it to [MATLAB's path](https://uk.mathworks.com/help/matlab/ref/pathtool.html). Requires MATLAB >= 2019a. Also needs [npy-matlab](https://github.com/kwikteam/npy-matlab), [brewermap](https://github.com/DrosteEffect/BrewerMap), and [prettify-matlab](https://github.com/Julie-Fabre/prettify_matlab).

**Both:** You need the [Allen CCF atlas files](https://github.com/cortex-lab/allenCCF) accessible locally.

```{toctree}
:maxdepth: 1
:caption: Contents

getting_started
find_and_fetch
plot_2d
plot_3d
plot_injections
compare_regions
thresholding
striatum_analysis
examples
api
contributing
```
