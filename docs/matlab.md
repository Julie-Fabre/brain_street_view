# MATLAB usage

Brain Street View was originally written in MATLAB. All MATLAB functions live in the `+bsv/` namespace package and are called as `bsv.functionName()`.

## Requirements

- MATLAB >= 2019a
- [allenCCF](https://github.com/cortex-lab/allenCCF) — Allen Atlas files
- [npy-matlab](https://github.com/kwikteam/npy-matlab) — NPY file I/O
- [brewermap](https://github.com/DrosteEffect/BrewerMap) — colormap generation
- [prettify-matlab](https://github.com/Julie-Fabre/prettify_matlab) — plot styling

## Installation

Clone the repository and add it (plus the dependencies above) to [MATLAB's path](https://uk.mathworks.com/help/matlab/ref/pathtool.html).

## Example workflow

```matlab
% 1. Find experiments
experimentIDs = bsv.findConnectivityExperiments({'VISp', 'VISl'});

% 2. Fetch data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, ...
    saveLocation, fileName, 'injectionIntensity', false, '', allenAtlasPath);

% 3. Plot 2D projections
bsv.plotConnectivity(experimentImgs, allenAtlasPath, 'CP', 10, 15, ...
    'coronal', true, 2, 'global', [], 'injectionIntensity')

% 4. 3D visualization
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, 'CP', ...
    'plotPatch', true)

% 5. Threshold
bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, 'CP', ...
    'percentile', 95)

% 6. Subregion analysis
bsv.analyzeCPSubregions(experimentImgs, allenAtlasPath, 'CP', ...
    'exportCSV', true)
```

See `+bsv/example.m` for the full workflow and `gettingStarted.mlx` for an interactive guide.

## Function reference

| Function | Description |
|---|---|
| `bsv.findConnectivityExperiments` | Query Allen API for experiments by region and mouse line |
| `bsv.fetchConnectivityData` | Download and cache projection density maps |
| `bsv.plotConnectivity` | 2D coronal/sagittal slice visualization |
| `bsv.plotConnectivity3D` | 3D isosurface rendering |
| `bsv.plotConnectivityMultiRegion` | Multi-region comparison plots |
| `bsv.thresholdConnectivity` | Apply signal thresholding |
| `bsv.analyzeCPSubregions` | Anatomical subdivision analysis with CSV export |
