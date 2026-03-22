import matplotlib.pyplot as plt

from .plot_connectivity import plot_connectivity


def plot_connectivity_multi_region(experiment_data, allen_atlas_path, output_regions,
                                   number_of_chunks, number_of_pixels, plane,
                                   region_only, smoothing, color_limits, color,
                                   normalization_info='unknown', input_regions=None,
                                   region_groups=None, experiment_region_info=None,
                                   normalize_by_group=False, custom_slices=None,
                                   slice_averaging=0):
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
