import os
import numpy as np
import bsv  # Assuming you'll create a Python module for bsv functions

# Define paths and parameters
save_location = '/home/julie/Dropbox/Data/AllenQueries'
allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF'
file_name = ''  # leave empty to recompute each time, or enter text to save and reload

# Input information about the regions and experiments you want to plot
input_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
mouse_line = ''  # leave empty to include all
primary_injection = True

# Define your plotting parameters
subtract_other_hemisphere = False
normalization_method = 'injectionIntensity'  # can be 'none' or 'injectionIntensity'
number_of_slices = 10
number_of_pixels = 15
output_regions = ['CP']
color = np.array([
    [0.543, 0, 0],
    [0, 0.746, 1],
    [0.180, 0.543, 0.340],
    [1, 0.547, 0]
])
plane = 'coronal'
smoothing = 2
color_limits = 'global'
region_only = True

# Get all Allen connectivity experiments of interest
experiment_ids = bsv.find_connectivity_experiments(input_regions, mouse_line, primary_injection)

# Fetch/load experiment data
experiment_imgs, injection_summary = bsv.fetch_connectivity_data(
    experiment_ids, save_location, file_name, normalization_method, subtract_other_hemisphere
)

# Plot projection data (in 2D)
bsv.plot_connectivity(
    experiment_imgs, allen_atlas_path, output_regions[0], number_of_slices, 
    number_of_pixels, plane, region_only, smoothing, color_limits, color
)

# Plot injections
# 2D, region by region
for input_region in input_regions:
    bsv.plot_connectivity(
        experiment_imgs, allen_atlas_path, input_region, number_of_slices, 
        number_of_pixels, plane, region_only, smoothing, color_limits, color
    )

# 3D, all regions
plot_patch = True  # if True plots a full volume; if False, plots a grid
bsv.plot_connectivity_3d(injection_summary, allen_atlas_path, output_regions[0], color, plot_patch)