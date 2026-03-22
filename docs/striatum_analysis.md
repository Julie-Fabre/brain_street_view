# Striatum subregion analysis

Break down projection density by anatomical subdivisions of the caudate putamen (CP) and nucleus accumbens (NAc) using the Allen v2 atlas, with optional CSV export.

```{image} ../images/docs/striatum_subregions.png
:width: 100%
```
*Mean fluorescence intensity per striatum subregion for projections from all visual cortical areas to CP and NAc. Subregions are defined by the Allen v2 atlas.*

Requires the Allen v2 atlas files (`annotation_volume_v2_20um_by_index.npy` and `UnifiedAtlas_Label_ontology_v2.csv`) in a separate directory.

**Python:**
```python
# First get the projection matrix from plot_connectivity
proj_matrix, proj_coords = bsv.plot_connectivity(
    experiment_data=experiment_imgs,
    allen_atlas_path='/path/to/allenCCF',
    output_region='CP',
    number_of_chunks=10,
    number_of_pixels=15,
    plane='coronal',
    region_only=True,
    smoothing=2,
    color_limits='global',
    color=None,
    normalization_info='injectionIntensity')

# Analyze subregions
subregion_results, global_results = bsv.analyze_cp_subregions(
    projection_matrix_array=proj_matrix,
    projection_matrix_coordinates_ara=proj_coords,
    allen_atlas_path_v2='/path/to/allenCCF_v2')

# Export to CSV
subregion_results, global_results = bsv.analyze_cp_subregions(
    projection_matrix_array=proj_matrix,
    projection_matrix_coordinates_ara=proj_coords,
    allen_atlas_path_v2='/path/to/allenCCF_v2',
    save_csv_path='/path/to/output.csv')

# Or do plot + analysis in one step
proj_matrix, proj_coords, sub_results, glob_results = \
    bsv.plot_connectivity_with_subregion_analysis(
        experiment_data=experiment_imgs,
        allen_atlas_path='/path/to/allenCCF',
        allen_atlas_path_v2='/path/to/allenCCF_v2',
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

**MATLAB:**
```matlab
[projMatrix, projCoords] = bsv.plotConnectivity(experimentImgs, allenAtlasPath, ...
    'CP', 10, 15, 'coronal', true, 2, 'global', [], 'injectionIntensity');

[subResults, globalResults] = bsv.analyzeCPSubregions(projMatrix, projCoords, ...
    allenAtlasPath_v2, [], inputRegions, [], csvPath);
```
