%% Neuro street view example script 

% local save directory 
saveLocation = '/home/julie/Dropbox/Data/AllenQueries';
allenAtlasPath =  '/home/julie/Dropbox/Atlas/allenCCF';
fileName = ''; % leave empty to recompute each time, or enter text (e.g. fileName = 'Visual_projections') to save and reload

% experiment to load information
inputRegions = {'VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor'}; % region(s), use Allen Atlas abbreviation conventions
    % or if you already have experiment IDs , skip step 1 below and input
    % them directly in step 2
mouseLine = ''; % leave empty to include all. use allen mouse line ids. 0 = wild-type. 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

% experiment loading parameters 
subtractOtherHemisphere = false;
normalizationMethod = 'injectionIntensity'; % can be 'none' or 'injectionIntensity'

% plotting parameters
numberOfSlices = 10;
numberOfPixels = 15;
outputRegions = {'CP'};
color = [0.543,0, 0; ...
    0, 0.746, 1;...
    0.180,0.543,0.340;...
    1,0.547,0]; % - not implemented yet - outline color for each region in RGB. leave empty to use defaults 
plane = 'coronal'; % - not implemented yet - coronal or sagital
smoothing = 2; % - not implemented yet - none or a number (of pixels)
colorLimits = 'global'; % - not implemented yet - global, per slice or two numbers  
regionOnly = true; % - not implemented yet - 

%% 1. Get allen connectivity experiments of interest 
experimentIDs = bsv.findConnectivityExperiments(inputRegions, mouseLine, primaryInjection);

%% 2. Fetch/load experiment data 
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, normalizationMethod, subtractOtherHemisphere, '', allenAtlasPath, false);

%% 3. Plot projection data (2D) 
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationMethod)

%% 4. Plot injections
%% a. 2D, region by region
for iInputRegion = 1:size(inputRegions,2)
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationMethod)
end

%% b. 2D, all regions - QQ TO DO
% bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

%% c. 3D, all regions
plotPatch = true; % if true plots a full volume; if false, plots a grid
bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, outputRegions(1), color, plotPatch)

