
%% Enhanced Neuro Street View Pipeline Script

% local save directory
username = char(java.lang.System.getProperty('user.name'));

saveLocation = ['/home/' username '/Dropbox/Data/AllenQueries'];
allenAtlasPath = ['/home/' username '/Dropbox/Atlas/allenCCF'];
fileName = ''; % leave empty to recompute each time, or enter text (e.g. fileName = 'Visual_projections') to save and reload

% experiment to load information
mouseLines = {'0', '177838435', '265180449', '177839022', '177836119'}; % multiple mouse lines to analyze
% wild type, Rbp4-Cre_KL100, Tlx3-Cre_PL56, Cux2-IRES-Cre, Emx1-IRES-Cre 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

% experiment loading parameters
subtractOtherHemisphere = false;
dataFetchNormalization = 'injectionVolume'; % can be 'none' or 'injectionVolume'

% plotting parameters
numberOfSlices = 10;
numberOfPixels = 15;
color = [0.543, 0, 0; ...
    0, 0.746, 1; ...
    0.180, 0.543, 0.340; ...
    1, 0.547, 0]; % outline color for each region in RGB
plane = 'coronal'; % coronal or sagital
smoothing = 2; % smoothing parameter
colorLimits = 'global'; % global, per slice or two numbers
outputRegions = {'CP', 'GPe', 'SNr'};
regionOnly = true;

%% 🔧 ENHANCED PROCESSING PARAMETERS
% These parameters improve cross-region comparison and analysis

% Threshold parameters
thresholdValue = 95; % For percentile: 95th percentile; for relative: 0.95; for zscore: 2.0; for absolute: 0.2
thresholdMethod = 'percentile'; % 'percentile', 'zscore', 'relative', 'absolute'

% Normalization parameters for enhanced comparison
analysisNormalization = 'region'; % 'region', 'zscore', 'robust', 'none'

fprintf('\n🚀 PIPELINE CONFIGURATION\n');
fprintf('Mouse lines: %s\n', strjoin(mouseLines, ', '));
fprintf('Threshold method: %s (value: %g)\n', thresholdMethod, thresholdValue);
fprintf('Analysis normalization: %s\n', analysisNormalization);
fprintf('Output regions: %s\n', strjoin(outputRegions, ', '));
fprintf('══════════════════════════════════\n');

%% ~~ VIS to CP ~~
inputRegions = {'VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

% Get allen connectivity experiments of interest
% Collect experiments from all mouse lines
allExperimentIDs = [];
for iMouseLine = 1:length(mouseLines)
    experimentIDs_thisLine = bsv.findConnectivityExperiments(inputRegions, mouseLines{iMouseLine}, primaryInjection);
    allExperimentIDs = [allExperimentIDs, experimentIDs_thisLine];
    fprintf('Mouse line %s: %d experiments\n', mouseLines{iMouseLine}, length(experimentIDs_thisLine));
end
experimentIDs = allExperimentIDs;
fprintf('Total experiments: %d\n', length(experimentIDs));

% Fetch/load experiment data
groupingMethod = '';
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, ...
    subtractOtherHemisphere, groupingMethod, allenAtlasPath, false);

% Plot projection data (2D) with enhanced visualization
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
    regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

% Define a "visual zone" with enhanced thresholding
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
    regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization, dataFetchNormalization);

save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/cp.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/cp_ara.mat', 'projectionMatrixCoordinates_ARA')


% Plot projection data (2D)
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, plane, ...
    regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

% Define a "visual zone"
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, plane, ...
    regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

% % Plot injections
% % a. 2D, region by region
% for iInputRegion = 1:size(inputRegions,2)
%     bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)
% end

% b. 2D, all regions - QQ TO DO
% bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

% % c. 3D, all regions
% plotPatch = true; % if true plots a full volume; if false, plots a grid
% bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, outputRegions(1), color, plotPatch)

%% ~~ DMS to GPe/SNr ~~
%manually selected expIDs
experimentIDs = [124059700, 159223001, 293366741, 293366035, 127762867, 159941339, 293366741, 301180385];
%160540013, ... % d2 301180385, ... % efr3a 112458831, ... % bit more
%ventral? % 292959343 - chat - beautiful!

%Fetch/load experiment data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, subtractOtherHemisphere, ...
    groupingMethod, allenAtlasPath, false);

%Plot projection data (2D)
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

% Define a "visual zone"
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization, dataFetchNormalization);


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe_ara.mat', 'projectionMatrixCoordinates_ARA')

[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization, dataFetchNormalization);


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr_ara.mat', 'projectionMatrixCoordinates_ARA')


[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, 99, 'percentile', analysisNormalization, dataFetchNormalization); % High threshold (99th percentile)


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr_thresh1.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr_ara_thresh1.mat', 'projectionMatrixCoordinates_ARA')

[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, 99, 'percentile', analysisNormalization, dataFetchNormalization); % High threshold (99th percentile)


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe_thresh1.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe_ara_thresh1.mat', 'projectionMatrixCoordinates_ARA')

% %Plot injections
% % a. 2D, region by region
% for iInputRegion = 1:size(inputRegions,2)
%     bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iInputRegion), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)
% end

% b. 2D, all regions - QQ TO DO
% bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

% % c. 3D, all regions
% plotPatch = true; % if true plots a full volume; if false, plots a grid
% bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, outputRegions(2), color, plotPatch)
% bsv.plotConnectivity3D(injectionSummary, allenAtlasPath, outputRegions(3), color, plotPatch)

%% ~~ DMS to GPe/SNr, by medial/lateral~~
%manually selected expIDs
experimentIDs = [124059700, 159223001, 293366741, 293366035, 127762867, 159941339, 293366741, 301180385];
%160540013, ... % d2 301180385, ... % efr3a 112458831, ... % bit more
%ventral? % 292959343 - chat - beautiful!

%Fetch/load experiment data
groupingMethod = 'AP';
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, subtractOtherHemisphere, ...
    groupingMethod, allenAtlasPath, false);

%Plot projection data (2D)
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

%%  ACA projections

inputRegions = {'ACA'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

% Get allen connectivity experiments of interest
% Collect experiments from all mouse lines
allExperimentIDs = [];
for iMouseLine = 1:length(mouseLines)
    experimentIDs_thisLine = bsv.findConnectivityExperiments(inputRegions, mouseLines{iMouseLine}, primaryInjection);
    allExperimentIDs = [allExperimentIDs, experimentIDs_thisLine];
    fprintf('Mouse line %s: %d experiments\n', mouseLines{iMouseLine}, length(experimentIDs_thisLine));
end
experimentIDs = allExperimentIDs;
fprintf('Total experiments: %d\n', length(experimentIDs));

% Fetch/load experiment data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, subtractOtherHemisphere, '', allenAtlasPath, false);

% Plot projection data (2D)
[projectionMatrix, projectionMatrix_coordinates] = bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), ...
    numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color);

% Define a "ACA zone"
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
    regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization);

save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/cp_aca.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/cp_ara_aca.mat', 'projectionMatrixCoordinates_ARA')


experimentIDs = [287994474, 159329308, 100142580, 127140981];
%160540013, ... % d2 301180385, ... % efr3a 112458831, ... % bit more
%ventral? % 292959343 - chat - beautiful!
groupingMethod = 'NaN';
%Fetch/load experiment data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, subtractOtherHemisphere, ...
    groupingMethod, allenAtlasPath, false);

%Plot projection data (2D)
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, dataFetchNormalization);

% Define a "visual zone"
[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(2), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization, dataFetchNormalization);


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe_aca.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/gpe_ara_aca.mat', 'projectionMatrixCoordinates_ARA')

[projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
    bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, outputRegions(3), numberOfSlices, numberOfPixels, ...
    plane, regionOnly, smoothing, colorLimits, color, thresholdValue, thresholdMethod, analysisNormalization, dataFetchNormalization);


save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr_aca.mat', 'projectionMatrix_array')
save('/home/julie/Dropbox/MATLAB/onPaths/brain_street_view/MATLAB/+bsv/bggRegions/visual/snr_ara_aca.mat', 'projectionMatrixCoordinates_ARA')

%% aud

inputRegions = {'AUDp'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

% Get allen connectivity experiments of interest
% Collect experiments from all mouse lines
allExperimentIDs = [];
for iMouseLine = 1:length(mouseLines)
    experimentIDs_thisLine = bsv.findConnectivityExperiments(inputRegions, mouseLines{iMouseLine}, primaryInjection);
    allExperimentIDs = [allExperimentIDs, experimentIDs_thisLine];
    fprintf('Mouse line %s: %d experiments\n', mouseLines{iMouseLine}, length(experimentIDs_thisLine));
end
experimentIDs = allExperimentIDs;
fprintf('Total experiments: %d\n', length(experimentIDs));

% Fetch/load experiment data
[experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, subtractOtherHemisphere, '', allenAtlasPath, false);

% Plot projection data (2D)
[projectionMatrix, projectionMatrix_coordinates] = bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), ...
    numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color);

%% different VIS to CP

inputRegions = {'VISp', 'VISl', 'VISal', 'VISam', 'VISpl', 'VISpm', 'VISli', 'VISpor', 'LP', 'LGd', 'LGv','PF' 'CL', 'IAM'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

for iInputRegions = 1:size(inputRegions, 2)
    
    % Get allen connectivity experiments of interest
    experimentIDs = bsv.findConnectivityExperiments(inputRegions(iInputRegions), mouseLine, primaryInjection);

    % Fetch/load experiment data
    groupingMethod = '';
    [experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, ...
        subtractOtherHemisphere, groupingMethod, allenAtlasPath, false);


    % Plot projection data (2D)
    [projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
        bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
        regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
end

%% and aud 

inputRegions = {'AUDp', 'AUDd', 'AUDv', 'MGd', 'MGv', 'MGm'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

for iInputRegions = 5:size(inputRegions, 2)
    
    % Get allen connectivity experiments of interest
    experimentIDs = bsv.findConnectivityExperiments(inputRegions(iInputRegions), mouseLine, primaryInjection);

    % Fetch/load experiment data
    groupingMethod = '';
    [experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, ...
        subtractOtherHemisphere, groupingMethod, allenAtlasPath, false);


    % Plot projection data (2D)
    [projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
        bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
        regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
end
% other senses 
inputRegions = {'OLF', 'VISC'}; % region(s), use Allen Atlas abbreviation conventions
% or if you already have experiment IDs , skip step 1 below and input
% them directly in step 2

for iInputRegions = 1:size(inputRegions, 2)
    
    % Get allen connectivity experiments of interest
    experimentIDs = bsv.findConnectivityExperiments(inputRegions(iInputRegions), mouseLine, primaryInjection);

    % Fetch/load experiment data
    groupingMethod = '';
    [experimentImgs, injectionSummary] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, dataFetchNormalization, ...
        subtractOtherHemisphere, groupingMethod, allenAtlasPath, false);


    % Plot projection data (2D)
    [projectionMatrix_array, projectionMatrixCoordinates_ARA] = ...
        bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions(1), numberOfSlices, numberOfPixels, plane, ...
        regionOnly, smoothing, colorLimits, color, dataFetchNormalization);
end

%% 🔬 ENHANCED THRESHOLDING COMPARISON DEMO
% Demonstrate different thresholding methods on the same data for comparison

fprintf('\n🔬 THRESHOLDING METHOD COMPARISON\n');
fprintf('Comparing different methods on last loaded data...\n');

% Use the last loaded data for comparison
if exist('experimentImgs', 'var') && exist('experimentIDs', 'var')
    testRegion = outputRegions{1}; % CP
    
    % Method 1: 95th Percentile
    fprintf('\n1️⃣ 95th Percentile Method:\n');
    [~, ~] = bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, testRegion, numberOfSlices, numberOfPixels, ...
        plane, regionOnly, smoothing, colorLimits, color, 95, 'percentile', 'region', dataFetchNormalization);
    
    % Method 2: 2-sigma Z-score  
    fprintf('\n2️⃣ Z-score Method (2σ):\n');
    [~, ~] = bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, testRegion, numberOfSlices, numberOfPixels, ...
        plane, regionOnly, smoothing, colorLimits, color, 2, 'zscore', 'zscore', dataFetchNormalization);
    
    % Method 3: 80% Relative
    fprintf('\n3️⃣ Relative Method (80%% of max):\n');
    [~, ~] = bsv.thresholdConnectivity(experimentImgs, allenAtlasPath, testRegion, numberOfSlices, numberOfPixels, ...
        plane, regionOnly, smoothing, colorLimits, color, 0.8, 'relative', 'none', dataFetchNormalization);
    
    fprintf('\n✅ Comparison complete! Check the generated figures to see differences.\n');
    fprintf('💡 TIP: Use ''percentile'' method for consistent cross-region comparisons.\n');
else
    fprintf('No experiment data loaded for comparison demo.\n');
end

fprintf('\n🎯 PIPELINE COMPLETE\n');
fprintf('📋 SUMMARY OF ENHANCEMENTS:\n');
fprintf('  🧬 Multiple mouse lines: %s\n', strjoin(mouseLines, ', '));
fprintf('  📊 Enhanced thresholding: %s method\n', thresholdMethod);
fprintf('  🎨 Professional visualization with colormaps\n');
fprintf('  📈 Cross-region statistical comparison\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

