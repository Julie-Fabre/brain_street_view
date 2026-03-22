import matplotlib.pyplot as plt

from .plot_connectivity import plot_connectivity


def plot_connectivity_multi_region(experiment_data, allen_atlas_path, output_regions,
                                   number_of_chunks, number_of_pixels, plane,
                                   region_only, smoothing, color_limits, color,
                                   normalization_info='unknown', input_regions=None,
                                   region_groups=None, experiment_region_info=None,
                                   normalize_by_group=False, custom_slices=None,
                                   slice_averaging=0):
    """Plot projections to multiple target regions for side-by-side comparison.

    Parameters
    ----------
    experiment_data : numpy.ndarray
        Projection density array from :func:`fetch_connectivity_data`.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    output_regions : list of str
        Target region acronyms (e.g. ``['CP', 'ACB']``).
    number_of_chunks : int
        Number of slices per region.
    number_of_pixels : int
        Pixel resolution per slice panel.
    plane : str
        ``'coronal'`` or ``'sagittal'``.
    region_only : bool
        Mask display to each target region boundary.
    smoothing : float
        Gaussian smoothing sigma in pixels.
    color_limits : str or list
        Colour scale specification.
    color : list or None
        RGB colour(s).
    normalization_info : str, optional
        Normalization label.
    input_regions : list of str, optional
        Source region acronyms for grouped display.
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
    proj_matrix_array : dict
        ``{region: proj_array}`` for each output region.
    proj_coords_array : dict
        ``{region: proj_coords}`` for each output region.
    """
    if isinstance(output_regions, str):
        output_regions = [output_regions]

    n_output_regions = len(output_regions)
    proj_matrix_array = {}
    proj_coords_array = {}

    print(f'Plotting {n_output_regions} output regions in combined view...')

    for i, region in enumerate(output_regions):
        print(f'Processing region {i + 1}/{n_output_regions}: {region}')

        proj_mat, proj_coords = plot_connectivity(
            experiment_data, allen_atlas_path, region,
            number_of_chunks, number_of_pixels, plane,
            region_only, smoothing, color_limits, color,
            normalization_info, input_regions, region_groups,
            experiment_region_info, normalize_by_group,
            custom_slices, slice_averaging)

        proj_matrix_array[region] = proj_mat
        proj_coords_array[region] = proj_coords

    return proj_matrix_array, proj_coords_array
