# Striatum subregion analysis

Break down projection density by anatomical subdivisions of the caudate putamen (CP) and nucleus accumbens (NAc) using the Allen v2 atlas. This generates per-slice and global statistics, with optional CSV export.

## Python

```python
import bsv

# First, plot projections and get the projection matrix
proj_matrix, proj_coords = bsv.plot_connectivity(
    experiment_imgs, '/path/to/allenCCF', 'CP',
    10, 15, 'coronal', True, 2, 'global', color,
    'injectionIntensity')

# Analyze CP and NAc subregions
subregion_results, global_results = bsv.analyze_cp_subregions(
    proj_matrix, proj_coords, '/path/to/allenCCF_v2')

# Export results to CSV
subregion_results, global_results = bsv.analyze_cp_subregions(
    proj_matrix, proj_coords, '/path/to/allenCCF_v2',
    input_regions=['VISp', 'VISl'],
    save_csv_path='/path/to/output/subregion_analysis.csv')

# Combined plot + analysis in one step
proj_matrix, proj_coords, subregion_results, global_results = \
    bsv.plot_connectivity_with_subregion_analysis(
        experiment_imgs, '/path/to/allenCCF', '/path/to/allenCCF_v2',
        'CP', 10, 15, 'coronal', True, 2, 'global', color,
        'injectionIntensity')
```

## MATLAB

```matlab
% Plot and get projection matrix
[projectionMatrix, projCoords] = bsv.plotConnectivity(experimentImgs, ...
    allenAtlasPath, 'CP', 10, 15, 'coronal', true, 2, 'global', color, ...
    'injectionIntensity');

% Analyze subregions
[subregionResults, globalResults] = bsv.analyzeCPSubregions(projectionMatrix, ...
    projCoords, allenAtlasPath_v2);

% Export to CSV
[subregionResults, globalResults] = bsv.analyzeCPSubregions(projectionMatrix, ...
    projCoords, allenAtlasPath_v2, [], inputRegions, [], csvPath);

% Combined plot + analysis
[projMatrix, projCoords, subResults, globalResults] = ...
    bsv.plotConnectivityWithSubregionAnalysis(experimentImgs, allenAtlasPath, ...
    allenAtlasPath_v2, 'CP', 10, 15, 'coronal', true, 2, 'global', color, ...
    'injectionIntensity');
```

## Requirements

Striatum subregion analysis requires the **Allen v2 atlas** files:
- `annotation_volume_v2_20um_by_index.npy`
- `UnifiedAtlas_Label_ontology_v2.csv`

These are separate from the standard Allen CCF files and should be placed in their own directory.
