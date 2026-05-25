import numpy as np
import pandas as pd
import os
import urllib.request

# Auto-download URLs for atlas files.
# The annotation volume (~2.4 GB) is hosted on figshare:
#   https://figshare.com/articles/dataset/Modified_Allen_CCF_2017_for_cortex-lab_allenCCF/25365829
# The structure tree CSV is the small lookup table from cortex-lab/allenCCF.
# These are downloaded automatically the first time load_atlas() is called if the
# files are not already present at allen_atlas_path.
_AUTO_DOWNLOAD_URLS = {
    'annotation_volume_10um_by_index.npy': 'https://ndownloader.figshare.com/files/44925493',
    'structure_tree_safe_2017.csv': (
        'https://raw.githubusercontent.com/cortex-lab/allenCCF/master/structure_tree_safe_2017.csv'
    ),
}


def _download_file(url, dest_path, filename):
    """Download a file with a simple progress indicator."""
    print(f"Downloading {filename} ...")

    def _progress(block_num, block_size, total_size):
        if total_size > 0:
            downloaded = min(block_num * block_size, total_size)
            pct = downloaded * 100 // total_size
            mb_done = downloaded / 1e6
            mb_total = total_size / 1e6
            print(f"\r  {pct:3d}%  {mb_done:.0f} / {mb_total:.0f} MB", end='', flush=True)

    urllib.request.urlretrieve(url, dest_path, _progress)
    print(f"\r  Done. Saved to {dest_path}          ")


def get_atlas_files(atlas_type='allen', atlas_resolution=10):
    """Return annotation volume and structure tree filenames for an atlas.

    Parameters
    ----------
    atlas_type : str, optional
        ``'allen'`` or a custom atlas name.
    atlas_resolution : int, optional
        Resolution in micrometres (10 or 20 for Allen).

    Returns
    -------
    tuple of str
        ``(annotation_filename, structure_tree_filename)``.
    """
    if atlas_type.lower() == 'allen':
        if atlas_resolution == 10:
            return 'annotation_volume_10um_by_index.npy', 'structure_tree_safe_2017.csv'
        elif atlas_resolution == 20:
            return 'annotation_volume_v2_20um_by_index.npy', 'UnifiedAtlas_Label_ontology_v2.csv'
        else:
            raise ValueError(f'Unsupported Allen atlas resolution: {atlas_resolution}. Supported: 10, 20 um.')
    else:
        return f'{atlas_type}_annotation_{atlas_resolution}um.npy', f'{atlas_type}_structure_tree.csv'


def load_atlas(allen_atlas_path, atlas_type='allen', atlas_resolution=10):
    """Load the annotation volume and structure tree for an atlas.

    Parameters
    ----------
    allen_atlas_path : str
        Path to the atlas directory.
    atlas_type : str, optional
        ``'allen'`` or a custom atlas name.
    atlas_resolution : int, optional
        Resolution in micrometres (10 or 20 for Allen).

    Returns
    -------
    av : numpy.ndarray
        Annotation volume (AP x DV x ML).
    st : pandas.DataFrame
        Structure tree.
    """
    annotation_file, structure_file = get_atlas_files(atlas_type, atlas_resolution)
    allen_atlas_path = str(allen_atlas_path)
    os.makedirs(allen_atlas_path, exist_ok=True)

    for filename in (annotation_file, structure_file):
        filepath = os.path.join(allen_atlas_path, filename)
        if not os.path.exists(filepath):
            if filename in _AUTO_DOWNLOAD_URLS:
                _download_file(_AUTO_DOWNLOAD_URLS[filename], filepath, filename)
            else:
                raise FileNotFoundError(
                    f"Atlas file not found and no auto-download URL is configured: {filepath}\n"
                    "Please download it manually and place it in: {allen_atlas_path}"
                )

    av = np.load(os.path.join(allen_atlas_path, annotation_file))
    st = load_structure_tree(os.path.join(allen_atlas_path, structure_file))
    return av, st


def load_structure_tree(filepath):
    """Load a structure tree CSV as a DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.read_csv(filepath)


def find_structure_indices(st, region_name):
    """Find structure indices matching a region name (prefix match like MATLAB version).
    Returns 1-based indices matching the annotation volume 'by_index' convention."""
    matches = st.index[st['acronym'].str.contains(region_name, na=False)].tolist()
    # Filter for exact prefix match, return 1-based (MATLAB find() convention)
    keep = []
    for idx in matches:
        acronym = st.loc[idx, 'acronym']
        if acronym[:len(region_name)] == region_name:
            keep.append(idx + 1)  # 1-based for annotation volume
    return keep


def get_structure_color(st, structure_idx):
    """Get RGB color [0-1] from hex triplet in structure tree.
    structure_idx is 1-based (from find_structure_indices), convert to 0-based for pandas."""
    hex_str = str(st.loc[structure_idx - 1, 'color_hex_triplet']).strip()
    # Pad to 6 chars if needed
    hex_str = hex_str.zfill(6)
    r = int(hex_str[0:2], 16) / 255.0
    g = int(hex_str[2:4], 16) / 255.0
    b = int(hex_str[4:6], 16) / 255.0
    return np.array([r, g, b])


def load_projection_info():
    """Load allenAtlasProjection_info.csv from docs/ relative to this package.
    Normalizes column names from hyphens to underscores (MATLAB VariableNamingRule='modify')."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(pkg_dir)
    csv_path = os.path.join(parent_dir, 'docs', 'allenAtlasProjection_info.csv')
    df = pd.read_csv(csv_path)
    df.columns = [c.replace('-', '_') for c in df.columns]
    return df
