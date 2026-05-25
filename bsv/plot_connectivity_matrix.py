"""Connectivity matrix heatmap visualization."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

from .atlas_utils import load_atlas, find_structure_indices


def plot_connectivity_matrix(
    projection_data,
    allen_atlas_path,
    target_regions,
    metric='mean',
    atlas_type='allen',
    atlas_resolution=10,
    normalize_rows=False,
    normalize_cols=False,
    cmap='viridis',
    annotate=True,
    figsize=None,
    title=None,
):
    """Create an N x M connectivity matrix heatmap.

    Build a matrix showing connectivity strength from N source regions to
    M target regions and display it as a heatmap.

    Parameters
    ----------
    projection_data : dict
        Dictionary mapping source region acronyms to tuples returned by
        ``fetch_connectivity_data``: ``(combined_projection, combined_info,
        individual_projections, experiment_info)``. Only the first element
        (the projection array) is used.
    allen_atlas_path : str
        Path to the directory containing atlas annotation volume and
        structure tree files.
    target_regions : list of str
        List of target region acronyms (e.g., ``['CP', 'ACB', 'SNr']``).
    metric : str, optional
        Connectivity metric to compute:

        - ``'mean'`` (default): Mean projection density in target region.
        - ``'max'``: Maximum projection density in target region.
        - ``'volume'``: Total projection volume (sum of density * voxel volume).
    atlas_type : str, optional
        Atlas type (default ``'allen'``).
    atlas_resolution : int, optional
        Atlas resolution in micrometres (default 10).
    normalize_rows : bool, optional
        If True, normalize each row (source) to [0, 1] by dividing by row max.
    normalize_cols : bool, optional
        If True, normalize each column (target) to [0, 1] by dividing by col max.
    cmap : str, optional
        Matplotlib colormap name (default ``'viridis'``).
    annotate : bool, optional
        If True, display numeric values in each cell (default True).
    figsize : tuple of float, optional
        Figure size ``(width, height)``. If None, auto-calculated from matrix size.
    title : str, optional
        Plot title. If None, no title is displayed.

    Returns
    -------
    connectivity_matrix : numpy.ndarray
        2D array of shape ``(n_sources, n_targets)`` with connectivity values.
    fig : matplotlib.figure.Figure
        The matplotlib Figure object.

    Examples
    --------
    >>> from bsv import fetch_connectivity_data, find_connectivity_experiments
    >>> from bsv import plot_connectivity_matrix
    >>> source_regions = ['VISam', 'VISp']
    >>> projection_data = {}
    >>> for src in source_regions:
    ...     exp_ids = find_connectivity_experiments(injection_structure=src)
    ...     data = fetch_connectivity_data(
    ...         experiment_ids=exp_ids,
    ...         save_location='./cache',
    ...         file_name=f'{src}_proj'
    ...     )
    ...     projection_data[src] = data
    >>> target_regions = ['CP', 'ACB']
    >>> matrix, fig = plot_connectivity_matrix(
    ...     projection_data=projection_data,
    ...     allen_atlas_path='./atlas',
    ...     target_regions=target_regions,
    ...     metric='mean'
    ... )
    """
    if metric not in ('mean', 'max', 'volume'):
        raise ValueError(f"metric must be 'mean', 'max', or 'volume', got '{metric}'")

    # Load atlas
    av, st = load_atlas(allen_atlas_path, atlas_type=atlas_type, atlas_resolution=atlas_resolution)

    # Projection grid is 100 um resolution: shape (132, 80, 114)
    # Atlas is typically 10 um: need to downsample atlas to match
    proj_resolution = 100  # um
    scale_factor = atlas_resolution / proj_resolution

    # Downsample annotation volume to projection grid resolution
    if scale_factor != 1.0:
        av_downsampled = zoom(av, scale_factor, order=0)  # nearest-neighbor
    else:
        av_downsampled = av

    # Build masks for each target region
    target_masks = {}
    for region in target_regions:
        indices = find_structure_indices(st, region)
        if not indices:
            print(f"Warning: No structures found matching '{region}'")
            target_masks[region] = np.zeros(av_downsampled.shape, dtype=bool)
        else:
            mask = np.isin(av_downsampled, indices)
            target_masks[region] = mask

    # Get source regions in consistent order
    source_regions = list(projection_data.keys())
    n_sources = len(source_regions)
    n_targets = len(target_regions)

    # Voxel volume in mm^3 (100 um = 0.1 mm per side)
    voxel_volume = (proj_resolution / 1000.0) ** 3

    # Build connectivity matrix
    connectivity_matrix = np.zeros((n_sources, n_targets))

    for i, src in enumerate(source_regions):
        data = projection_data[src]
        # Handle both tuple returns and direct array
        if isinstance(data, tuple):
            proj = data[0]  # combined_projection
        else:
            proj = data

        # Average across groups if multiple groups exist
        if proj.ndim == 4:
            proj = np.mean(proj, axis=3)

        for j, tgt in enumerate(target_regions):
            mask = target_masks[tgt]

            # Handle shape mismatch between projection and mask
            if proj.shape != mask.shape:
                # Resize mask to match projection
                zoom_factors = [p / m for p, m in zip(proj.shape, mask.shape)]
                mask_resized = zoom(mask.astype(float), zoom_factors, order=0) > 0.5
            else:
                mask_resized = mask

            masked_values = proj[mask_resized]

            if len(masked_values) == 0 or np.all(np.isnan(masked_values)):
                connectivity_matrix[i, j] = 0.0
            elif metric == 'mean':
                connectivity_matrix[i, j] = np.nanmean(masked_values)
            elif metric == 'max':
                connectivity_matrix[i, j] = np.nanmax(masked_values)
            elif metric == 'volume':
                connectivity_matrix[i, j] = np.nansum(masked_values) * voxel_volume

    # Apply normalization
    if normalize_rows:
        row_max = connectivity_matrix.max(axis=1, keepdims=True)
        row_max[row_max == 0] = 1  # Avoid division by zero
        connectivity_matrix = connectivity_matrix / row_max

    if normalize_cols:
        col_max = connectivity_matrix.max(axis=0, keepdims=True)
        col_max[col_max == 0] = 1
        connectivity_matrix = connectivity_matrix / col_max

    # Create figure
    if figsize is None:
        figsize = (max(6, n_targets * 1.2), max(4, n_sources * 0.8))

    fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    im = ax.imshow(connectivity_matrix, cmap=cmap, aspect='auto')

    # Set ticks and labels
    ax.set_xticks(np.arange(n_targets))
    ax.set_yticks(np.arange(n_sources))
    ax.set_xticklabels(target_regions)
    ax.set_yticklabels(source_regions)

    # Rotate x labels for readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    # Add colorbar
    metric_labels = {
        'mean': 'Mean projection density',
        'max': 'Max projection density',
        'volume': 'Total projection volume (mm³)',
    }
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(metric_labels[metric])

    # Annotate cells with values
    if annotate:
        # Determine text color based on background
        norm_matrix = (connectivity_matrix - connectivity_matrix.min())
        if connectivity_matrix.max() > connectivity_matrix.min():
            norm_matrix = norm_matrix / (connectivity_matrix.max() - connectivity_matrix.min())

        for i in range(n_sources):
            for j in range(n_targets):
                value = connectivity_matrix[i, j]
                # Use white text on dark cells, black on light
                text_color = 'white' if norm_matrix[i, j] < 0.5 else 'black'
                # Format based on magnitude
                if metric == 'volume':
                    text = f'{value:.2e}' if value < 0.01 else f'{value:.3f}'
                else:
                    text = f'{value:.2e}' if value < 0.001 else f'{value:.4f}'
                ax.text(j, i, text, ha='center', va='center', color=text_color, fontsize=8)

    # Labels and title
    ax.set_xlabel('Target Region')
    ax.set_ylabel('Source Region')
    if title:
        ax.set_title(title)

    plt.tight_layout()

    return connectivity_matrix, fig
