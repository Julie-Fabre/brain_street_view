"""Generate striatum subregion analysis figure for docs."""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'axes.linewidth': 1.5,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import bsv

save_location = '/home/julie/Dropbox/Data/AllenQueries'
allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF'
allen_atlas_path_v2 = '/home/julie/Dropbox/Atlas/allenCCF_v2'
img_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images', 'docs')

print('Finding experiments...')
input_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
exp_ids = bsv.find_connectivity_experiments(input_regions)

print('Fetching data...')
exp_imgs, inj_summary, _, _ = bsv.fetch_connectivity_data(
    exp_ids, save_location, '', 'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path)

print('Plotting connectivity...')
proj_array, proj_coords = bsv.plot_connectivity(
    exp_imgs, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity')
plt.close('all')

print('Running subregion analysis...')
sub_results, glob_results = bsv.analyze_cp_subregions(
    proj_array, proj_coords, allen_atlas_path_v2)

out_path = os.path.join(img_dir, 'striatum_subregions.png')
plt.savefig(out_path, facecolor='white', edgecolor='none')
plt.close('all')
print(f'Saved: {out_path}')
