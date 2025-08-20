%% Example: Using custom atlas with brain_street_view
% This example demonstrates how to use different atlas types and resolutions

%% Example 1: Using Allen atlas with 20um resolution
allenAtlasPath = '/path/to/your/atlas';
experimentIDs = [123456, 789012]; % Your experiment IDs
saveLocation = './connectivity_data';
fileName = 'my_connectivity_data';

% Fetch data with 20um resolution Allen atlas
[combinedProjection, combinedInjectionInfo, ~, experimentRegionInfo] = ...
    bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    'injectionVolume', false, 'NaN', allenAtlasPath, false, [], [], false, false, ...
    'allen', 20); % Specify atlas type and resolution

% Plot connectivity with 20um resolution
outputRegion = 'CP'; % Caudate putamen
[projectionMatrix, coordinates] = bsv.plotConnectivity(combinedProjection, ...
    allenAtlasPath, outputRegion, 4, 100, 'coronal', true, 2, [], 'viridis', ...
    'injectionVolume', [], [], experimentRegionInfo, false, [], 0, ...
    'allen', 20); % Specify atlas type and resolution

%% Example 2: Using a custom atlas
% For custom atlases, place your atlas files with naming convention:
% - [atlasType]_annotation_[resolution]um.npy
% - [atlasType]_structure_tree.csv

customAtlasPath = '/path/to/custom/atlas';
atlasType = 'myatlas'; % Your custom atlas name
atlasResolution = 25; % Your atlas resolution in microns

% Ensure your atlas files are named:
% - myatlas_annotation_25um.npy
% - myatlas_structure_tree.csv

% Fetch data with custom atlas
[combinedProjection, combinedInjectionInfo, ~, experimentRegionInfo] = ...
    bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    'injectionVolume', false, 'NaN', customAtlasPath, false, [], [], false, false, ...
    atlasType, atlasResolution);

% Plot with custom atlas
[projectionMatrix, coordinates] = bsv.plotConnectivity(combinedProjection, ...
    customAtlasPath, outputRegion, 4, 100, 'coronal', true, 2, [], 'viridis', ...
    'injectionVolume', [], [], experimentRegionInfo, false, [], 0, ...
    atlasType, atlasResolution);

%% Example 3: 3D visualization with custom resolution
injectionSummary = combinedInjectionInfo; % Your injection data
regionToPlot = 'CP';
color = [1 0 0]; % Red
plotPatch = true;

% 3D plot with custom atlas parameters
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, regionToPlot, ...
    color, plotPatch, 'allen', 20); % Using 20um Allen atlas

%% Notes:
% 1. The atlas type and resolution parameters are optional and default to:
%    - atlasType: 'allen'
%    - atlasResolution: 10
%
% 2. For Allen atlas, supported resolutions are:
%    - 10um: Uses annotation_volume_10um_by_index.npy and structure_tree_safe_2017.csv
%    - 20um: Uses annotation_volume_v2_20um_by_index.npy and UnifiedAtlas_Label_ontology_v2.csv
%
% 3. For custom atlases, ensure your files follow the naming convention:
%    - Annotation: [atlasType]_annotation_[resolution]um.npy
%    - Structure tree: [atlasType]_structure_tree.csv
%
% 4. The saved data files will include atlas information in their names:
%    - [fileName]_[normalization]_sub[subtract]_[atlasType]_[resolution]um.mat