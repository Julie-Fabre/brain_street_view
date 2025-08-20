%% Example: Analyzing CP Subregions and Exporting to CSV
% This example demonstrates how to use analyzeCPSubregions with CSV export

%% Step 1: Run plotConnectivity to get projection data
allenAtlasPath = '/path/to/your/allenCCF';
allenAtlasPath_v2 = '/path/to/your/allenCCF_v2'; % v2 atlas for subregion analysis

% Your connectivity data from fetchConnectivityData
combinedProjection = ...; % Your projection data
outputRegion = 'CP'; % Analyze caudate putamen

% Run plotConnectivity
[projectionMatrix, projectionCoordinates] = bsv.plotConnectivity(...
    combinedProjection, allenAtlasPath, outputRegion, ...
    4, 100, 'coronal', true, 2, [], 'viridis');

%% Step 2: Analyze CP subregions WITHOUT CSV export
[subregionResults, globalResults] = bsv.analyzeCPSubregions(...
    projectionMatrix, projectionCoordinates, allenAtlasPath_v2);

%% Step 3: Analyze CP subregions WITH CSV export
csvSavePath = './CP_subregion_analysis.csv';
[subregionResults, globalResults] = bsv.analyzeCPSubregions(...
    projectionMatrix, projectionCoordinates, allenAtlasPath_v2, ...
    [], [], [], csvSavePath);

%% Step 4: Analyze with specific parameters and CSV export
% Specify slices, input regions, and groups
outputSlices = [1, 2, 3]; % Analyze specific slices
inputRegions = {'MOp', 'MOs', 'SS'}; % Motor and sensory regions
regionGroups = [1, 1, 2]; % Group motor regions together

csvSavePath = './CP_subregion_analysis_grouped.csv';
[subregionResults, globalResults] = bsv.analyzeCPSubregions(...
    projectionMatrix, projectionCoordinates, allenAtlasPath_v2, ...
    outputSlices, inputRegions, regionGroups, csvSavePath);

%% What the CSV contains:
% The exported CSV file will include:
% 1. SubregionID - Unique ID for each CP subregion
% 2. SubregionName - Full name of the subregion
% 3. SubregionAcronym - Abbreviated name
% 4. GlobalMeanIntensity - Average intensity across all groups/slices
% 5. TotalVoxelCount - Total number of voxels analyzed
% 6. MeanIntensity_[GroupName] - Mean intensity for each group
% 7. VoxelCount_[GroupName] - Voxel count for each group

% Additionally, a second CSV file with '_nonzero' suffix will be created
% containing only subregions with non-zero mean intensities.

%% Example of loading and working with the CSV
% Load the CSV file
results = readtable(csvSavePath);

% Find top 5 subregions by global mean intensity
sortedResults = sortrows(results, 'GlobalMeanIntensity', 'descend');
top5Subregions = sortedResults(1:min(5, height(sortedResults)), :);

fprintf('Top 5 CP subregions by mean intensity:\n');
for i = 1:height(top5Subregions)
    fprintf('%d. %s: %.4f\n', i, ...
        top5Subregions.SubregionAcronym{i}, ...
        top5Subregions.GlobalMeanIntensity(i));
end

%% Notes:
% - The CSV export is optional - leave saveCsvPath empty to skip export
% - Two CSV files are created: full results and non-zero results only
% - Column names are automatically generated based on input region names
% - Results are sorted by global mean intensity (highest first)
% - The CSV can be easily imported into Excel, R, or Python for further analysis