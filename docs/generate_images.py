"""Generate documentation images with publication-quality styling."""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

# Pretty defaults: sans-serif, larger text, thicker lines for publication
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 20,
    'axes.titlesize': 24,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.linewidth': 3.5,
    'axes.linewidth': 2.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import bsv

# Paths - update these for your system
save_location = '/tmp/AllenQueries'
allen_atlas_path = '/home/jf5479/Dropbox/Atlas/allenCCF'
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

# ── 2. AP-grouped: VISam → CP split by injection AP location ──
print('Generating: VISam → CP grouped by AP location...')
exp_imgs_ap, _, _, exp_region_info_ap = bsv.fetch_connectivity_data(
    exp_ids, save_location, '', 'injectionIntensity', False,
    allen_atlas_path=allen_atlas_path,
    grouping_method='AP')

bsv.plot_connectivity(
    exp_imgs_ap, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity')
save_current('plot_ap_grouped_VISam_CP.png')

# ── 3. Plot injection sites ──
print('Generating: VISam injection sites...')
bsv.plot_connectivity(
    exp_imgs, allen_atlas_path, 'VISam',
    10, 15, 'coronal', True, 2, 'global', None, 'injectionIntensity')
save_current('plot_injections_VISam.png')

# ── 4. Thresholded connectivity ──
print('Generating: thresholded VISam → CP...')
bsv.threshold_connectivity(
    exp_imgs, allen_atlas_path, 'CP',
    10, 15, 'coronal', True, 2, 'global', None,
    threshold=90, threshold_method='percentile',
    normalization_method='none', data_fetch_normalization='injectionIntensity')
save_current('threshold_VISam_CP.png')

# ── 5. Multi-region grouped: all visual areas → CP ──
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

# ── 6. 3D visualization (static frame) ──
print('Generating: 3D VISam → CP...')
bsv.plot_connectivity_3d(inj_summary, allen_atlas_path, 'CP',
                          plot_patch=True, animate=False)
save_current('plot_3d_VISam_CP.png')

# ── 7. Connectivity matrix heatmap ──
print('Generating: connectivity matrix heatmap...')
source_regions = ['VISp', 'VISl', 'VISam', 'VISpm']
target_regions = ['CP', 'ACB', 'SNr']

projection_data = {}
for src in source_regions:
    src_ids = bsv.find_connectivity_experiments([src])
    src_data = bsv.fetch_connectivity_data(
        src_ids, save_location, '', 'injectionIntensity', False,
        allen_atlas_path=allen_atlas_path)
    projection_data[src] = src_data

matrix, fig = bsv.plot_connectivity_matrix(
    projection_data=projection_data,
    allen_atlas_path=allen_atlas_path,
    target_regions=target_regions,
    metric='mean',
    cmap='viridis',
    annotate=True,
    title='Visual → Subcortical Connectivity')
save_current('plot_connectivity_matrix.png')

print('\nDone! All images in:', img_dir)
