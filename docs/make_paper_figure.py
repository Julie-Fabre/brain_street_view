"""Create composite figure for JOSS paper."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import string

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 14,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

img_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images', 'docs')

img_2d = mpimg.imread(os.path.join(img_dir, 'plot_2d_VISam_CP.png'))
img_grouped = mpimg.imread(os.path.join(img_dir, 'plot_grouped_VIS_CP.png'))
img_3d = mpimg.imread(os.path.join(img_dir, 'plot_3d_VISam_CP.png'))

fig = plt.figure(figsize=(16, 12), facecolor='white')

# A: 2D projections (top left, wide)
ax_a = fig.add_axes([0.02, 0.68, 0.96, 0.30])
ax_a.imshow(img_2d)
ax_a.axis('off')
ax_a.text(-0.01, 1.05, 'A', transform=ax_a.transAxes,
          fontsize=22, fontweight='bold', va='top', ha='left')

# B: Grouped regions (bottom left)
ax_b = fig.add_axes([0.02, 0.02, 0.60, 0.62])
ax_b.imshow(img_grouped)
ax_b.axis('off')
ax_b.text(-0.01, 1.05, 'B', transform=ax_b.transAxes,
          fontsize=22, fontweight='bold', va='top', ha='left')

# C: 3D (bottom right)
ax_c = fig.add_axes([0.62, 0.05, 0.37, 0.56])
ax_c.imshow(img_3d)
ax_c.axis('off')
ax_c.text(-0.03, 1.05, 'C', transform=ax_c.transAxes,
          fontsize=22, fontweight='bold', va='top', ha='left')

out_path = os.path.join(img_dir, 'fig_paper.png')
plt.savefig(out_path, facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out_path}')
