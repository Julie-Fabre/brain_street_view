import os
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display, clear_output

from .atlas_utils import load_atlas, find_structure_indices, load_projection_info
from .fetch_connectivity_images import fetch_connectivity_images, PROJECTION_GRID_SIZE
from .fetch_connectivity_summary import fetch_connectivity_summary, load_injection_summary


def plot_upstream_projectome(experiment_ids, source_regions, target_region, save_location,
                              allen_atlas_path, atlas_resolution=10, atlas_type='allen',
                              slice_thickness=2, density_percentile=99,
                              static_ap=None, save_path=None):
    """Interactive coronal slice viewer of upstream projections into a target region.

    Displays a Jupyter widget with a coronal AP slider showing:

    - Atlas structure boundaries rendered at full atlas resolution
    - Per-source-region coloured heatmaps of projection density in the target
      region voxels (intensity ∝ projection density, colour = source region)
    - Injection site scatter dots sized by injection volume

    Parameters
    ----------
    experiment_ids : list of int
        Experiment IDs, as returned by :func:`find_connectivity_experiments`
        with *target_regions* set.
    source_regions : list of str
        Source region acronyms to display (subset of upstream regions).
    target_region : str
        Target brain region acronym (e.g. ``'CP'``).
    save_location : str
        Directory where ``density.raw`` and injection summary files are cached.
    allen_atlas_path : str
        Path to the Allen CCF atlas directory.
    atlas_resolution : int, optional
        Atlas resolution in micrometres (default 10).
    atlas_type : str, optional
        Atlas type (default ``'allen'``).
    slice_thickness : int, optional
        Number of 100 µm grid slices averaged either side of the displayed AP
        index when computing the projection heatmap (default 2).
    density_percentile : float, optional
        Percentile of non-zero target-voxel density values used as the colour
        saturation ceiling (default 99). Lower values make dim signals brighter;
        raise toward 100 to avoid saturation on bright spots.
    static_ap : int or None, optional
        If given, render this single AP slice (100 µm grid index) as a static
        figure instead of launching the interactive widget — used to produce
        documentation figures headlessly.
    save_path : str or None, optional
        When *static_ap* is set, save the rendered figure to this path.
    """
    AP, DV, ML = PROJECTION_GRID_SIZE
    # Number of full-res atlas voxels per 100 µm projection grid voxel
    factor = 100 // atlas_resolution

    # ------------------------------------------------------------------ atlas
    print("Loading atlas...")
    av, st = load_atlas(allen_atlas_path, atlas_type, atlas_resolution)

    # Full-resolution atlas dimensions (DV x ML for display)
    DV_full = av.shape[1]
    ML_full = av.shape[2]

    target_idx = find_structure_indices(st, target_region)

    # Target mask at 100 µm grid resolution (for density computation)
    av_grid = av[::factor, ::factor, ::factor][:AP, :DV, :ML]
    target_mask = np.isin(av_grid, target_idx)   # (AP, DV, ML) bool
    del av_grid  # free memory — only needed for target_mask

    # ---------------------------------------------------------- experiment map
    projection_info = load_projection_info()
    exp_to_region = {}
    for exp_id in experiment_ids:
        row = projection_info[projection_info['id'] == exp_id]
        if len(row) > 0:
            abbrev = row.iloc[0]['structure_abbrev']
            # Match exact or parent prefix so 'ACAd' maps to source region 'ACA'
            for src in source_regions:
                if abbrev == src or abbrev.startswith(src):
                    exp_to_region[exp_id] = src
                    break

    valid_ids = [e for e in experiment_ids if e in exp_to_region]

    # ------------------------------------------------------- colours (tab10)
    cmap = plt.cm.get_cmap('tab10', max(len(source_regions), 1))
    region_colors = {reg: np.array(cmap(i)[:3]) for i, reg in enumerate(source_regions)}

    # ----------------------------------------- load density volumes & metadata
    print(f"Loading {len(valid_ids)} experiment volumes...")
    densities = {}   # exp_id → (AP, DV, ML) float32
    inj_ap_coord = {}   # in atlas voxel units (atlas_resolution µm each)
    inj_dv_coord = {}
    inj_ml_coord = {}
    inj_vol = {}

    for exp_id in valid_ids:
        exp_dir = os.path.join(save_location, str(exp_id))

        raw_path = os.path.join(exp_dir, 'density.raw')
        if not os.path.exists(raw_path):
            fetch_connectivity_images(exp_id, exp_dir)
        if os.path.exists(raw_path):
            vol = np.fromfile(raw_path, dtype='<f4').reshape(PROJECTION_GRID_SIZE, order='F')
            densities[exp_id] = vol

        row = projection_info[projection_info['id'] == exp_id]
        if len(row) > 0:
            r = row.iloc[0]
            try:
                coords_str = str(r['injection_coordinates'])
                coords = [float(x) for x in coords_str.strip('[]').split(',')]
                # CSV coordinates are in µm; convert to atlas voxel indices
                inj_ap_coord[exp_id] = coords[0] / atlas_resolution
                inj_dv_coord[exp_id] = coords[1] / atlas_resolution
                inj_ml_coord[exp_id] = coords[2] / atlas_resolution
            except Exception:
                inj_ap_coord[exp_id] = inj_dv_coord[exp_id] = inj_ml_coord[exp_id] = -1
            try:
                inj_vol[exp_id] = float(r['injection_volume'])
            except Exception:
                inj_vol[exp_id] = 0.1

    # -------------------------------- pre-compute mean density per source region
    region_density = {}
    for reg in source_regions:
        reg_ids = [e for e in valid_ids if exp_to_region.get(e) == reg and e in densities]
        if reg_ids:
            stack = np.stack([densities[e] for e in reg_ids], axis=0)
            region_density[reg] = np.mean(stack, axis=0)
        else:
            region_density[reg] = np.zeros(PROJECTION_GRID_SIZE, dtype=np.float32)

    # Normalisation ceiling: 95th percentile of non-zero target-voxel densities
    # across all source regions, so dim signals are still clearly visible.
    all_target_vals = np.concatenate([
        region_density[reg][target_mask]
        for reg in source_regions if target_mask.any()
    ]) if source_regions else np.array([])
    nonzero = all_target_vals[all_target_vals > 0]
    global_max = float(np.percentile(nonzero, density_percentile)) if nonzero.size else 1.0
    if global_max == 0:
        global_max = 1.0

    # ---------------------------------------------- dot-size scaling (20-200)
    vols_arr = np.array([inj_vol.get(e, 0.1) for e in valid_ids])
    v_min, v_max = vols_arr.min(), vols_arr.max()
    if v_max == v_min:
        v_max = v_min + 1e-9
    dot_sizes = 20 + 180 * (vols_arr - v_min) / (v_max - v_min)

    # ------------------------------- full-resolution atlas background helper
    def _atlas_bg(ap_idx):
        """(DV_full, ML_full) grayscale at full atlas resolution."""
        # Use the middle atlas slice of the 100 µm window
        ap_atlas = min(ap_idx * factor + factor // 2, av.shape[0] - 1)
        sl = av[ap_atlas]       # (DV_full, ML_full) structure indices

        brain = sl > 0
        target = np.isin(sl, target_idx)

        h = sl[:, :-1] != sl[:, 1:]
        v = sl[:-1, :] != sl[1:, :]
        bnd = np.zeros(sl.shape, dtype=bool)
        bnd[:, :-1] |= h;  bnd[:, 1:] |= h
        bnd[:-1, :] |= v;  bnd[1:, :] |= v
        bnd &= brain

        # White brain, light-gray target region, black boundaries + outside
        img = np.zeros(sl.shape, dtype=float)   # black outside brain
        img[brain] = 1.0                         # white brain tissue
        img[target] = 0.85                       # light gray target region
        img[bnd] = 0.0                           # black structure boundaries
        return img

    # ----------------------------------------------------------- draw function
    def _draw(ap_idx):
        fig, ax = plt.subplots(figsize=(9, 5))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')

        # Full-resolution atlas background (DV_full x ML_full)
        ax.imshow(_atlas_bg(ap_idx), cmap='gray', vmin=0, vmax=1,
                  origin='upper', aspect='equal', interpolation='nearest',
                  extent=[0, ML_full, DV_full, 0])

        # Coloured heatmaps — density at 100 µm, upsampled to full atlas res
        ap_lo = max(0, ap_idx - slice_thickness)
        ap_hi = min(AP - 1, ap_idx + slice_thickness)
        tm = target_mask[ap_idx]  # (DV, ML) at 100 µm

        sorted_regions = sorted(
            source_regions,
            key=lambda r: float(np.mean(region_density[r][ap_idx][tm])) if tm.any() else 0.0,
        )
        for reg in sorted_regions:
            color = region_colors[reg]
            density_slice = np.mean(region_density[reg][ap_lo:ap_hi + 1], axis=0)  # (DV, ML)

            alpha = np.clip(density_slice / global_max, 0, 1)
            alpha[~tm] = 0.0

            # Upsample to full atlas resolution then trim to exact atlas shape
            alpha_up = np.repeat(np.repeat(alpha, factor, axis=0), factor, axis=1)
            alpha_up = alpha_up[:DV_full, :ML_full]

            rgba = np.zeros((DV_full, ML_full, 4), dtype=float)
            rgba[:, :, :3] = color
            rgba[:, :, 3] = alpha_up

            ax.imshow(rgba, origin='upper', aspect='equal', interpolation='nearest',
                      extent=[0, ML_full, DV_full, 0])

        # Injection site scatter (coords in atlas voxel units)
        ap_center = ap_idx * factor + factor / 2.0
        ap_window = slice_thickness * factor + factor  # half-window in atlas voxels
        for i, exp_id in enumerate(valid_ids):
            ap_val = inj_ap_coord.get(exp_id, -1)
            if abs(ap_val - ap_center) > ap_window:
                continue
            reg = exp_to_region.get(exp_id)
            if reg is None:
                continue
            ml_val = inj_ml_coord.get(exp_id, -1)
            dv_val = inj_dv_coord.get(exp_id, -1)
            if ml_val < 0 or dv_val < 0:
                continue
            ax.scatter(ml_val, dv_val, s=dot_sizes[i],
                       c=[region_colors[reg]],
                       edgecolors='white', linewidths=0.5,
                       zorder=10, alpha=0.9)

        # Legend
        handles = [
            plt.scatter([], [], s=80, c=[region_colors[r]],
                        edgecolors='white', linewidths=0.5, label=r)
            for r in source_regions
        ]
        if handles:
            leg = ax.legend(handles=handles, loc='upper right',
                            framealpha=0.3, labelcolor='white', fontsize=9)
            leg.get_frame().set_facecolor('black')

        ax.set_xlim(0, ML_full)
        ax.set_ylim(DV_full, 0)
        ax.set_xticks([]);  ax.set_yticks([])
        ax.set_title(
            f'Upstream projections → {target_region}   |   AP {ap_idx}  ({ap_idx * 100} µm)',
            color='white', fontsize=11, pad=6)

        plt.tight_layout()
        plt.show()

    # ------------------------------------- static render (docs / headless use)
    if static_ap is not None:
        _draw(int(static_ap))
        if save_path:
            plt.savefig(save_path, facecolor='black', bbox_inches='tight')
            print(f'Saved: {save_path}')
        return

    # --------------------------------------------------------- widget assembly
    out = widgets.Output()

    slider = widgets.IntSlider(
        min=0, max=AP - 1, value=AP // 2,
        description='AP (100 µm)',
        continuous_update=False,
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='520px'),
    )

    def _on_change(change):
        with out:
            clear_output(wait=True)
            _draw(slider.value)

    slider.observe(_on_change, names='value')

    display(widgets.VBox([slider, out]))
    with out:
        _draw(slider.value)
