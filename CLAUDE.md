# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Brain Street View is a MATLAB toolbox for loading and plotting Allen Connectivity Data from the Allen Institute Brain Connectivity Atlas (Oh et al., Nature, 2014). The toolbox allows researchers to query, visualize, and analyze neural projection data from mouse brain connectivity experiments.

## Architecture

The codebase is organized into a MATLAB package structure under `+bsv/`:

### Core Functions
- `findConnectivityExperiments.m` - Queries Allen Brain Map API to find experiments matching region criteria
- `fetchConnectivityData.m` - Downloads and processes experiment data, handles normalization and grouping
- `fetchConnectivitySummary.m` - Retrieves injection site metadata and statistics
- `fetchConnectivityImages.m` - Downloads raw density data files from Allen Institute
- `plotConnectivity.m` - Creates 2D visualizations of connectivity data across brain slices
- `plotConnectivity3D.m` - Generates 3D visualizations of injection sites and projections
- `thresholdConnectivity.m` - Applies threshold filtering to connectivity data

### Data Flow
1. **Query** → `findConnectivityExperiments()` searches Allen API for experiments
2. **Fetch** → `fetchConnectivityData()` downloads raw data and metadata
3. **Process** → Data normalization, hemisphere subtraction, grouping
4. **Visualize** → `plotConnectivity()` or `plotConnectivity3D()` for analysis

### Key Dependencies
- allenCCF (Allen Common Coordinate Framework atlas files)
- npy-matlab (for reading NumPy files)
- brewermap (colormaps)
- prettify-matlab (plot formatting)

## Development Commands

This is a MATLAB-only project with no build system. Development workflow:

1. **Getting Started**: Run `gettingStarted.mlx` (MATLAB Live Script)
2. **Example Usage**: See `MATLAB/+bsv/example.m` for basic workflow
3. **Testing**: No automated test suite - verify functionality through example scripts

## Data Management

- Raw data cached locally in `saveLocation` directory structure
- First run downloads images (slow), subsequent runs load cached data (fast)
- Data organized by experiment ID in subdirectories
- Supports normalization by injection intensity/volume

## API Integration

The toolbox integrates with Allen Brain Map REST APIs:
- Experiment queries: `http://api.brain-map.org/api/v2/data/query.json`
- Data downloads: `http://api.brain-map.org/grid_data/download/`
- Metadata: `http://connectivity.brain-map.org/api/v2/data/ProjectionStructureUnionize/`

## Visualization Parameters

Key plotting parameters in `plotConnectivity()`:
- `numberOfSlices` - Number of coronal/sagittal sections
- `numberOfPixels` - Resolution of projection binning
- `plane` - 'coronal' or 'sagital' orientation
- `regionOnly` - Limit display to specific regions
- `normalizationMethod` - 'injectionIntensity', 'injectionVolume', or 'none'

## Special Notes

- Spikes should be plotted below "presence ratio" (user instruction)
- Data coordinates: AP (anterior-posterior), DV (dorsal-ventral), ML (medial-lateral)
- Hemisphere handling: Can subtract contralateral projections for specificity analysis
- Region specification uses Allen Atlas abbreviations (e.g., 'VISp', 'CP', 'SNr')