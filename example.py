"""Brain Street View example script - Python port of +bsv/example.m"""
import bsv

# Local paths
save_location = '/home/julie/Dropbox/Data/AllenQueries'
allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF'
file_name = ''  # leave empty to recompute each time

# Experiment parameters
input_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
mouse_line = ''
primary_injection = True

# Loading parameters
subtract_other_hemisphere = False
normalization_method = 'injectionIntensity'

# Plotting parameters
number_of_slices = 10
number_of_pixels = 15
output_regions = ['CP']
color = [[0.543, 0, 0], [0, 0.746, 1], [0.180, 0.543, 0.340], [1, 0.547, 0]]
plane = 'coronal'
smoothing = 2
color_limits = 'global'
region_only = True

# 1. Find experiments
experiment_ids = bsv.find_connectivity_experiments(input_regions, mouse_line, primary_injection)

# 2. Fetch/load experiment data
experiment_imgs, injection_summary, _, _ = bsv.fetch_connectivity_data(
    experiment_ids, save_location, file_name,
    normalization_method, subtract_other_hemisphere,
    allen_atlas_path=allen_atlas_path)

# 3. Plot projection data (2D)
bsv.plot_connectivity(experiment_imgs, allen_atlas_path, output_regions[0],
                       number_of_slices, number_of_pixels, plane,
                       region_only, smoothing, color_limits, color,
                       normalization_method)

# 4. Plot injections - region by region
for region in input_regions:
    bsv.plot_connectivity(experiment_imgs, allen_atlas_path, region,
                           number_of_slices, number_of_pixels, plane,
                           region_only, smoothing, color_limits, color,
                           normalization_method)

# 5. 3D plot
bsv.plot_connectivity_3d(injection_summary, allen_atlas_path, output_regions[0],
                          color, plot_patch=True)
