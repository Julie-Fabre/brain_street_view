%% Neuro street view example script 

% local save directory 
saveLocation = '/home/julie/Dropbox/Data/AllenQueries';

% experiment to load information
regions = {'VISp', 'VISl', 'VISal', 'VISpm', 'VISam'}; % region(s), use Allen Atlas abbreviation conventions
    % can also be a single or vector of experiment IDs e.g. experiments = [183618139, 112952510, 177459319];
mouseLine = ''; % leave empty to include all 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

% experiment loading parameters 
normalizationMethod = 'injectionIntensity'; % can be 'none', 'injectionVolume', 'injectionIntensity'
subtractOtherHemisphere = true;

% plotting parameters

%% Get allen connectivity experiments of interest 
experimentIDs = nsv_findConnectivityExperiments(regions, mouseLine, primaryInjection);

%% Fetch/load experiment data 
experimentData = nsv_fetchConnectivityData(experimentIDs, saveLocation, normalizationMethod, subtractOtherHemisphere);

%% Plot injections (2D)


%% Plot projection data (2D) 

