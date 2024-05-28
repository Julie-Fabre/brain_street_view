%% Neuro street view example script 

% local save directory 
saveLocation = '/home/julie/Dropbox/Data/AllenQueries';
allenAtlasPath =  '/home/julie/Dropbox/Atlas/allenCCF';
fileName = ''; % leave empty to recompute each time, or enter text (e.g. fileName = 'Visual_projections') to save and reload

% experiment to load information
inputRegions = {'VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor'}; % region(s), use Allen Atlas abbreviation conventions
    % or if you already have experiment IDs , skip step 1 below and input
    % them directly in step 2
mouseLine = ''; % leave empty to include all 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

% experiment loading parameters 
subtractOtherHemisphere = true;
normalizationMethod = 'injectionIntensity'; % can be 'none' or 'injectionIntensity'

% plotting parameters
numberOfSlices = 10;
numberOfPixels = 15;
color = [0.543,0, 0; ...
    0, 0.746, 1;...
    0.180,0.543,0.340;...
    1,0.547,0]; % - not implemented yet - outline color for each region in RGB. leave empty to use defaults 
plane = 'coronal'; % - not implemented yet - coronal or sagital
smoothing = 2; % - not implemented yet - none or a number (of pixels)
colorLimits = 'global'; % - not implemented yet - global, per slice or two numbers  
outputRegions = {'CP', 'GPe', 'SNr'};
regionOnly = true; % - not implemented yet - 

%% 1. Get allen connectivity experiments of interest 
experimentIDs = nsv_findConnectivityExperiments(inputRegions, mouseLine, primaryInjection);

%% 2. Fetch/load experiment data 
experimentData = nsv_fetchConnectivityData(experimentIDs, saveLocation, fileName, normalizationMethod, subtractOtherHemisphere);

%% 3. Plot injections
%% a. 2D, region by region
for iInputRegion = 1:size(inputRegions,2)
    nsv_plotConnectivity(experimentData, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)
end

%% b. 2D, all regions
nsv_plotConnectivity(experimentData, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

%% b. 3D, all regions
nsv_plotConnectivity(experimentData, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)


%% 4. Plot projection data (2D) 
for iInputRegion = 1:size(outputRegions,2)
    nsv_plotConnectivity(experimentData, allenAtlasPath, outputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)
end