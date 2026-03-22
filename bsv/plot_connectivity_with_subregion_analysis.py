from .plot_connectivity import plot_connectivity
from .analyze_cp_subregions import analyze_cp_subregions


def plot_connectivity_with_subregion_analysis(experiment_data, allen_atlas_path,
                                               allen_atlas_path_v2, output_region,
                                               number_of_chunks, number_of_pixels, plane,
                                               region_only, smoothing, color_limits, color,
                                               normalization_info='unknown', input_regions=None,
                                               region_groups=None, experiment_region_info=None,
                                               normalize_by_group=False, custom_slices=None,
                                               slice_averaging=0):
    """Plot projections and run CP subregion analysis in one step.

    Combines :func:`plot_connectivity` and :func:`analyze_cp_subregions`.

    Parameters
    ----------
    experiment_data : numpy.ndarray
        Projection density array from :func:`fetch_connectivity_data`.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    allen_atlas_path_v2 : str
        Path to the Allen v2 atlas directory.
    output_region : str
        Target region acronym (e.g. ``'CP'``).
    number_of_chunks : int
        Number of slices.
    number_of_pixels : int
        Pixel resolution per slice panel.
    plane : str
        ``'coronal'`` or ``'sagittal'``.
    region_only : bool
        Mask to the target region boundary.
    smoothing : float
        Gaussian smoothing sigma in pixels.
    color_limits : str or list
        Colour scale specification.
    color : list or None
        RGB colour(s).
    normalization_info : str, optional
        Normalization label.
    input_regions : list of str, optional
        Source region acronyms.
    region_groups : list of int, optional
        Group assignment per input region.
    experiment_region_info : dict, optional
        Per-experiment metadata.
    normalize_by_group : bool, optional
        Normalize each group independently.
    custom_slices : list of int, optional
        Specific slice indices.
    slice_averaging : int, optional
        Adjacent-slice averaging radius.

    Returns
    -------
    proj_array : numpy.ndarray
        Projection matrix.
    proj_coords : list
        Coordinate information.
    subregion_results : dict or None
        Subregion analysis results (None if region is not CP).
    global_results : dict or None
        Global subregion results (None if region is not CP).
    """
    print('=== CONNECTIVITY ANALYSIS WITH CP SUBREGION ANALYSIS ===')

    # Step 1: connectivity plotting
    print('Step 1: Running connectivity analysis...')
    proj_array, proj_coords = plot_connectivity(
        experiment_data, allen_atlas_path, output_region,
        number_of_chunks, number_of_pixels, plane,
        region_only, smoothing, color_limits, color,
        normalization_info, input_regions, region_groups,
        experiment_region_info, normalize_by_group,
        custom_slices, slice_averaging)

    # Step 2: subregion analysis
    subregion_results = None
    global_results = None

    region_str = output_region if isinstance(output_region, str) else output_region[0]
    if 'CP' in region_str.upper():
        print('\nStep 2: Analyzing CP subregions...')
        subregion_results, global_results = analyze_cp_subregions(
            proj_array, proj_coords, allen_atlas_path_v2,
            input_regions=input_regions, region_groups=region_groups)
    else:
        print(f'Warning: Output region is not CP - skipping subregion analysis')

    print('\n=== ANALYSIS COMPLETE ===')
    return proj_array, proj_coords, subregion_results, global_results
