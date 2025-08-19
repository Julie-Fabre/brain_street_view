function [projectionMatrix_array, projectionMatrixCoordinates_ARA, subregionResults, globalResults] = plotConnectivityWithSubregionAnalysis(experimentData, allenAtlasPath, allenAtlasPath_v2, outputRegion, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationInfo, inputRegions, regionGroups, experimentRegionInfo, normalizeByGroup, customSlices, sliceAveraging)
% plotConnectivityWithSubregionAnalysis - Plot connectivity and analyze CP subregions
%
% This function combines plotConnectivity with CP subregion analysis
%
% INPUTS: Same as plotConnectivity, plus:
%   allenAtlasPath_v2 - Path to v2 atlas files for subregion analysis
%
% OUTPUTS: Same as plotConnectivity, plus:
%   subregionResults - Per-slice subregion analysis
%   globalResults - Global (all slices) subregion analysis

% Handle optional parameters (same as plotConnectivity)
if nargin < 18 || isempty(sliceAveraging)
    sliceAveraging = 0;
end
if nargin < 17 || isempty(customSlices)
    customSlices = [];
end
if nargin < 16 || isempty(normalizeByGroup)
    normalizeByGroup = false;
end

fprintf('=== CONNECTIVITY ANALYSIS WITH CP SUBREGION ANALYSIS ===\n');

%% Step 1: Run normal connectivity plotting
fprintf('Step 1: Running connectivity analysis...\n');

[projectionMatrix_array, projectionMatrixCoordinates_ARA] = bsv.plotConnectivity(experimentData, allenAtlasPath, outputRegion, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationInfo, inputRegions, regionGroups, experimentRegionInfo, normalizeByGroup, customSlices, sliceAveraging);

%% Step 2: Analyze CP subregions
fprintf('\nStep 2: Analyzing CP subregions...\n');

if strcmpi(outputRegion, 'CP') || contains(upper(outputRegion), 'CP')
    [subregionResults, globalResults] = bsv.analyzeCPSubregions(projectionMatrix_array, projectionMatrixCoordinates_ARA, allenAtlasPath_v2, [], inputRegions, regionGroups);
    
    %% Step 3: Create summary table and save results
    fprintf('\nStep 3: Creating summary table...\n');
    
    % Create a summary table
    summaryTable = table();
    summaryTable.SubregionID = globalResults.subregion_ids;
    summaryTable.SubregionName = globalResults.subregion_names;
    summaryTable.SubregionAcronym = globalResults.subregion_acronyms;
    summaryTable.GlobalMeanIntensity = globalResults.mean_intensities;
    summaryTable.TotalVoxelCount = globalResults.total_voxel_counts;
    
    % Add per-slice statistics
    summaryTable.MeanAcrossSlices = mean(subregionResults.mean_intensities, 2, 'omitnan');
    summaryTable.StdAcrossSlices = std(subregionResults.mean_intensities, 0, 2, 'omitnan');
    summaryTable.SlicesWithData = sum(~isnan(subregionResults.mean_intensities), 2);
    
    % Remove subregions with no data
    hasData = ~isnan(summaryTable.GlobalMeanIntensity) | summaryTable.SlicesWithData > 0;
    summaryTable = summaryTable(hasData, :);
    
    % Sort by global mean intensity
    summaryTable = sortrows(summaryTable, 'GlobalMeanIntensity', 'descend', 'MissingPlacement', 'last');
    
    fprintf('\nCP Subregion Summary (top 10):\n');
    fprintf('%-30s %-15s %-15s %-10s\n', 'Subregion', 'Global Mean', 'Mean±Std', 'Slices');
    fprintf('%-30s %-15s %-15s %-10s\n', repmat('-', 1, 30), repmat('-', 1, 15), repmat('-', 1, 15), repmat('-', 1, 10));
    
    nToShow = min(10, height(summaryTable));
    for i = 1:nToShow
        if ~isnan(summaryTable.GlobalMeanIntensity(i))
            fprintf('%-30s %-15.4f %-15s %-10d\n', ...
                summaryTable.SubregionAcronym{i}, ...
                summaryTable.GlobalMeanIntensity(i), ...
                sprintf('%.4f±%.4f', summaryTable.MeanAcrossSlices(i), summaryTable.StdAcrossSlices(i)), ...
                summaryTable.SlicesWithData(i));
        end
    end
    
    % Optionally save the summary table
    % Uncomment the following lines if you want to save to CSV:
    % csvFilename = sprintf('CP_subregion_analysis_%s.csv', datestr(now, 'yyyymmdd_HHMMSS'));
    % writetable(summaryTable, csvFilename);
    % fprintf('\nSummary table saved to: %s\n', csvFilename);
    
else
    warning('Output region is not CP - skipping subregion analysis');
    subregionResults = [];
    globalResults = [];
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');

end