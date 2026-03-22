from .plot_connectivity import plot_connectivity


def plot_multi_region_injections(experiment_imgs, allen_atlas_path, input_regions,
                                 number_of_slices, number_of_pixels, plane,
                                 region_only, smoothing, color_limits, color,
                                 normalization_method, experiment_region_info=None):
    """Plot injection sites for each input region separately.

    Calls :func:`plot_connectivity` once per region in *input_regions*.

    Parameters
    ----------
    experiment_imgs : numpy.ndarray
        Projection density array from :func:`fetch_connectivity_data`.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    input_regions : list of str
        Source region acronyms.
    number_of_slices : int
        Number of slices per region.
    number_of_pixels : int
        Pixel resolution per slice panel.
    plane : str
        ``'coronal'`` or ``'sagittal'``.
    region_only : bool
        Mask to each region boundary.
    smoothing : float
        Gaussian smoothing sigma in pixels.
    color_limits : str or list
        Colour scale specification.
    color : list or None
        RGB colour(s).
    normalization_method : str
        Normalization label.
    experiment_region_info : dict, optional
        Per-experiment metadata.
    """
    for region in input_regions:
        print(f'Processing injection region: {region}')
        plot_connectivity(experiment_imgs, allen_atlas_path, region,
                          number_of_slices, number_of_pixels, plane,
                          region_only, smoothing, color_limits, color,
                          normalization_method)
