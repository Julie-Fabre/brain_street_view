%% Striatum subregion analysis example
% Analyze projection intensity per CP and NAc subregion using the Allen v2 atlas.

%% Parameters
saveLocation = '/home/julie/Dropbox/Data/AllenQueries';
allenAtlasPath = '/home/julie/Dropbox/Atlas/allenCCF';
allenAtlasPath_v2 = '/home/julie/Dropbox/Atlas/allenCCF_v2';
fileName = '';

inputRegions = {'VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor'};
mouseLine = '';
primaryInjection = true;

subtractOtherHemisphere = false;
normalizationMethod = 'injectionIntensity';

numberOfSlices = 10;
numberOfPixels = 15;
outputRegion = 'CP';
plane = 'coronal';
smoothing = 2;
colorLimits = 'global';
regionOnly = true;
color = [];

%% 1. Fetch data and plot projections to CP
experimentIDs = bsv.findConnectivityExperiments(inputRegions, mouseLine, primaryInjection);

[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    normalizationMethod, subtractOtherHemisphere, '', allenAtlasPath, false);

[projectionMatrix_array, projectionMatrixCoordinates_ARA] = bsv.plotConnectivity(experimentImgs, ...
    allenAtlasPath, outputRegion, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, ...
    colorLimits, color, normalizationMethod);

%% 2. CP and NAc subregion analysis
[subregionResults, globalResults] = bsv.analyzeCPSubregions(projectionMatrix_array, ...
    projectionMatrixCoordinates_ARA, allenAtlasPath_v2);

%% 3. Combined connectivity + subregion analysis
[projectionMatrix_array, projectionMatrixCoordinates_ARA, subregionResults, globalResults] = ...
    bsv.plotConnectivityWithSubregionAnalysis(experimentImgs, allenAtlasPath, allenAtlasPath_v2, ...
    outputRegion, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, ...
    color, normalizationMethod);

%% 4. Export results to CSV
csvPath = fullfile(saveLocation, 'cp_subregion_analysis.csv');
[subregionResults, globalResults] = bsv.analyzeCPSubregions(projectionMatrix_array, ...
    projectionMatrixCoordinates_ARA, allenAtlasPath_v2, [], inputRegions, [], csvPath);
fprintf('Results saved to: %s\n', csvPath);
