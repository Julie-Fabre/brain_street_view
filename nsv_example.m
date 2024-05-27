%% Neuro street view example script 

% local save directory 
saveLocation = '/home/julie/Dropbox/Data/AllenQueries';
allenAtlasPath =  '/home/julie/Dropbox/Atlas/allenCCF';
fileName = 'Vis_projections';

% experiment to load information
regions = {'VISp', 'VISl', 'VISal', 'VISpm', 'VISam'}; % region(s), use Allen Atlas abbreviation conventions
    % or if you already have experiment IDs , skip step 1 below and input
    % them directly in step 2
mouseLine = ''; % leave empty to include all 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

% experiment loading parameters 
normalizationMethod = 'injectionIntensity'; % can be 'none' or 'injectionIntensity'
subtractOtherHemisphere = true;

% plotting parameters
color 
numberOfSlices 
plane % coroncal or sagital 
smoothing 
colorLimits 

%% 1. Get allen connectivity experiments of interest 
experimentIDs = nsv_findConnectivityExperiments(regions, mouseLine, primaryInjection);

%% 2. Fetch/load experiment data 
experimentData = nsv_fetchConnectivityData(experimentIDs, saveLocation, fileName, normalizationMethod, subtractOtherHemisphere);

%% 3. Plot injections (2D)
nsv_plotConnectivity(experimentData, regions, numberOfSlices, )

%% 4. Plot projection data (2D) 
nsv_plotConnectivity(experimentData, brainRegionToPlot)
