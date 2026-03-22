from .plot_connectivity import plot_connectivity


def plot_multi_region_injections(experiment_imgs, allen_atlas_path, input_regions,
                                 number_of_slices, number_of_pixels, plane,
                                 region_only, smoothing, color_limits, color,
                                 normalization_method, experiment_region_info=None):
    """Plot injection sites for each input region by calling plot_connectivity per region."""
    for region in input_regions:
        print(f'Processing injection region: {region}')
        plot_connectivity(experiment_imgs, allen_atlas_path, region,
                          number_of_slices, number_of_pixels, plane,
                          region_only, smoothing, color_limits, color,
                          normalization_method)
