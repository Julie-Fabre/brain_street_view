# Threshold projections

Identify significant projection signals by applying thresholds to density maps. Several thresholding methods are available.

## Python

```python
import bsv

# Absolute threshold
bsv.threshold_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    threshold=0.5, threshold_method='absolute')

# Percentile threshold (e.g. top 5% of voxels)
bsv.threshold_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    threshold=95, threshold_method='percentile')

# Z-score threshold
bsv.threshold_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    threshold=2.0, threshold_method='zscore')

# Relative threshold (% of max value)
bsv.threshold_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    threshold=0.1, threshold_method='relative')
```

## MATLAB

```matlab
bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, 'CP', ...
    10, 15, 'coronal', true, 2, 'global', color, 95, 'percentile')
```

## Threshold methods

| Method | `threshold` value | Description |
|---|---|---|
| `'absolute'` | density value | Keep voxels above this raw value |
| `'percentile'` | 0--100 | Keep voxels above this percentile |
| `'zscore'` | z-score cutoff | Keep voxels with z-score above threshold |
| `'relative'` | 0--1 | Keep voxels above this fraction of the max |
