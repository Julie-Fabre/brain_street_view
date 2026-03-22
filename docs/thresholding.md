# Threshold projections

Identify significant projection signals using absolute, percentile, z-score, or relative thresholds.

```{image} ../images/docs/threshold_VISam_CP.png
:width: 100%
```
*VISam projections to CP thresholded at the 90th percentile. Red dots mark voxels above threshold.*

**Python:**
```python
# Percentile — keep top 10% of voxels
bsv.threshold_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    threshold=90, threshold_method='percentile')

# Other methods:
#   threshold_method='absolute'  — raw density cutoff
#   threshold_method='zscore'    — z-score cutoff
#   threshold_method='relative'  — fraction of max (0–1)
```

**MATLAB:**
```matlab
bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, 'CP', ...
    10, 15, 'coronal', true, 2, 'global', color, 90, 'percentile')
```
