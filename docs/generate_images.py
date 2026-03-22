"""Generate documentation images with publication-quality styling."""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

# Pretty defaults: Arial, larger text, slightly thicker lines
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

# Paths
save_location = '/home/julie/Dropbox/Data/AllenQueries'
allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF'
img_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images', 'docs')
os.makedirs(img_dir, exist_ok=True)


def save_current(name):
    path = os.path.join(img_dir, name)
    plt.savefig(path, facecolor='white', edgecolor='none')
    plt.close('all')
    print(f'  Saved: {path}')


# ── 1. Single region: VISam → CP (2D projections) ──
print('Generating: VISam → CP projections (2D)...')
exp_ids = bsv.find_connectivity_experiments(['VISam'])
exp_imgs, inj_summary, _, _ = bsv.fetch_connectivity_data(
    exp_ids, save_location, '', 'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path)

proj_array, proj_coords = bsv.plot_connectivity(
    exp_imgs, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity')
save_current('plot_2d_VISam_CP.png')

# ── 2. Plot injection sites ──
print('Generating: VISam injection sites...')
bsv.plot_connectivity(
    exp_imgs, allen_atlas_path, 'VISam',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity')
save_current('plot_injections_VISam.png')

# ── 3. Thresholded connectivity ──
print('Generating: thresholded VISam → CP...')
bsv.threshold_connectivity(
    exp_imgs, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None,
    threshold=90, threshold_method='percentile',
    normalization_method='none', data_fetch_normalization='injectionIntensity')
save_current('threshold_VISam_CP.png')

# ── 4. Multi-region grouped: all visual areas → CP ──
print('Generating: grouped visual areas → CP...')
grouped_regions = ['VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor']
region_groups = [1, 2, 2, 3, 3, 3, 3, 3]
grouped_ids = bsv.find_connectivity_experiments(grouped_regions)
exp_imgs_g, _, _, exp_region_info_g = bsv.fetch_connectivity_data(
    grouped_ids, save_location, '', 'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path,
    input_regions=grouped_regions, region_groups=region_groups)

bsv.plot_connectivity(
    exp_imgs_g, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity',
    input_regions=grouped_regions, region_groups=region_groups,
    experiment_region_info=exp_region_info_g)
save_current('plot_grouped_VIS_CP.png')

# ── 5. 3D visualization (static frame) ──
print('Generating: 3D VISam → CP...')
bsv.plot_connectivity_3d(inj_summary, allen_atlas_path, 'CP',
                          plot_patch=True, animate=False)
save_current('plot_3d_VISam_CP.png')

print('\nDone! All images in:', img_dir)
