function [subregionResults, globalResults] = analyzeCPSubregions(projectionMatrix_array, projectionMatrixCoordinates_ARA, allenAtlasPath_v2, outputSlices, inputRegions, regionGroups)
% analyzeCPSubregions - Calculate mean fluorescence intensity per CP subregion
%
% INPUTS:
%   projectionMatrix_array - 3D array from plotConnectivity [pixels x pixels x slices]  
%   projectionMatrixCoordinates_ARA - Coordinate information from plotConnectivity
%   allenAtlasPath_v2 - Path to v2 atlas files (e.g., '/home/user/Dropbox/Atlas/allenCCF_v2')
%   outputSlices - Optional: specific slices to analyze (default: all)
%   inputRegions - Optional: cell array of input region names for plot titles
%   regionGroups - Optional: numeric array indicating region grouping
%
% OUTPUTS:
%   subregionResults - Structure with per-slice results
%   globalResults - Structure with global (all slices combined) results

if nargin < 4 || isempty(outputSlices)
    outputSlices = 1:size(projectionMatrix_array, 3);
end
if nargin < 5 || isempty(inputRegions)
    inputRegions = [];
end
if nargin < 6 || isempty(regionGroups)
    regionGroups = [];
end

fprintf('\n=== CP SUBREGION ANALYSIS ===\n');
fprintf('Analyzing CP subregions using v2 atlas...\n');

%% Load v2 atlas files
fprintf('Loading v2 atlas files...\n');

% Load annotation volume (20um resolution)
av_v2_file = fullfile(allenAtlasPath_v2, 'annotation_volume_v2_20um_by_index.npy');
if ~exist(av_v2_file, 'file')
    error('Cannot find v2 annotation file: %s', av_v2_file);
end
av_v2 = readNPY(av_v2_file);

% Load ontology
ontology_v2_file = fullfile(allenAtlasPath_v2, 'UnifiedAtlas_Label_ontology_v2.csv');
if ~exist(ontology_v2_file, 'file')
    error('Cannot find v2 ontology file: %s', ontology_v2_file);
end
ontology_v2 = readtable(ontology_v2_file);

fprintf('Atlas dimensions: %s\n', mat2str(size(av_v2)));
fprintf('Ontology entries: %d\n', height(ontology_v2));

%% Find CP subregions in ontology
fprintf('Finding CP subregions in ontology...\n');

% Find all caudoputamen entries
cp_mask = contains(ontology_v2.name, 'Caudoputamen', 'IgnoreCase', true) | ...
          contains(ontology_v2.acronym, 'CP', 'IgnoreCase', true);

cp_subregions = ontology_v2(cp_mask, :);
fprintf('Found %d CP subregions:\n', height(cp_subregions));

for i = 1:height(cp_subregions)
    fprintf('  %d: %s (%s)\n', cp_subregions.id(i), cp_subregions.name{i}, cp_subregions.acronym{i});
end

%% Initialize results structures (will update after determining data structure)
nSubregions = height(cp_subregions);

globalResults = struct();
globalResults.subregion_ids = cp_subregions.id;
globalResults.subregion_names = cp_subregions.name;
globalResults.subregion_acronyms = cp_subregions.acronym;
globalResults.mean_intensities = NaN(nSubregions, 1);
globalResults.total_voxel_counts = zeros(nSubregions, 1);

%% Resolution conversion factors
% Atlas v2 is 20um resolution, original plotting was based on 10um
resolution_factor = 2; % 20um / 10um

fprintf('\nAnalyzing fluorescence per subregion...\n');

%% Global analysis (all slices combined)
fprintf('Computing global averages across all slices...\n');

global_intensity_sums = zeros(nSubregions, 1);
global_voxel_counts = zeros(nSubregions, 1);

%% Determine data structure
fprintf('Data dimensions: %s\n', mat2str(size(projectionMatrix_array)));
nCoords = length(projectionMatrixCoordinates_ARA);
fprintf('Number of coordinate sets: %d\n', nCoords);

if ndims(projectionMatrix_array) == 4
    % 4D case: [pixels, pixels, groups, slices]
    nActualSlices = size(projectionMatrix_array, 4);
    nGroups = size(projectionMatrix_array, 3);
    fprintf('Found %d groups and %d slices in 4D data\n', nGroups, nActualSlices);
elseif ndims(projectionMatrix_array) == 3
    % Check coordinate info to determine if 3rd dim is groups or slices
    thirdDim = size(projectionMatrix_array, 3);
    fprintf('Third dimension size: %d\n', thirdDim);
    
    if nCoords == 1 && thirdDim > 1
        % 1 coordinate set but multiple in 3rd dim = groups
        nActualSlices = 1;
        nGroups = thirdDim;
        fprintf('DETECTED: %d groups and %d slice in 3D data (groups in 3rd dimension)\n', nGroups, nActualSlices);
    elseif nCoords == thirdDim
        % Number of coords matches 3rd dim = slices
        nActualSlices = thirdDim;
        nGroups = 1;
        fprintf('DETECTED: %d slices and %d group in 3D data (slices in 3rd dimension)\n', nActualSlices, nGroups);
    elseif nCoords < thirdDim && nCoords > 1
        % Fewer coordinates than 3rd dim but more than 1 = groups sharing coordinates
        nActualSlices = 1;  % Treat as single slice with multiple groups
        nGroups = thirdDim;
        fprintf('DETECTED: %d groups and %d slice in 3D data (groups sharing coordinates)\n', nGroups, nActualSlices);
    else
        % Based on your feedback: assume groups in 3rd dimension when we have multiple items there
        if thirdDim > 1
            nActualSlices = 1;
            nGroups = thirdDim;
            fprintf('OVERRIDE: %d groups and %d slice in 3D data (override based on expected structure)\n', nGroups, nActualSlices);
        else
            nActualSlices = thirdDim;
            nGroups = 1;
            fprintf('DEFAULT: %d slices and %d group in 3D data\n', nActualSlices, nGroups);
        end
    end
else
    error('Unexpected data dimensions: %s', mat2str(size(projectionMatrix_array)));
end

fprintf('FINAL: nActualSlices=%d, nGroups=%d\n', nActualSlices, nGroups);

% Update outputSlices to match actual data
if nargin < 4 || isempty(outputSlices)
    outputSlices = 1:nActualSlices;
else
    outputSlices = outputSlices(outputSlices <= nActualSlices);
end
nSlices = length(outputSlices);

fprintf('Processing %d slices with %d groups: %s\n', nSlices, nGroups, mat2str(outputSlices));

%% Initialize results structures with correct dimensions
subregionResults = struct();
subregionResults.slice_numbers = outputSlices;
subregionResults.subregion_ids = cp_subregions.id;
subregionResults.subregion_names = cp_subregions.name;
subregionResults.subregion_acronyms = cp_subregions.acronym;
subregionResults.nGroups = nGroups;
% Results matrix: [subregions x (groups * slices)]
subregionResults.mean_intensities = NaN(nSubregions, nGroups * nSlices);
subregionResults.voxel_counts = zeros(nSubregions, nGroups * nSlices);
subregionResults.coordinates_ARA = projectionMatrixCoordinates_ARA;

%% Process each slice
if nGroups > 1 && nActualSlices == 1
    % Special case: Groups in 3rd dimension, only 1 slice
    fprintf('\nProcessing 1 slice with %d groups...\n', nGroups);
    
    % Use coordinates from slice 1
    coords = projectionMatrixCoordinates_ARA{1};
    fprintf('  Using coordinates from set 1\n');
elseif nGroups > nActualSlices
    % Special case: More groups than slices (groups sharing coordinates)
    fprintf('\nProcessing %d slices with %d groups (groups sharing coordinates)...\n', nActualSlices, nGroups);
    
    % For now, use coordinates from slice 1 and process all groups
    coords = projectionMatrixCoordinates_ARA{1};
    fprintf('  Using coordinates from set 1 for all groups\n');
    % Set nActualSlices to 1 to simplify processing
    nActualSlices = 1;
    nSlices = 1;
else
    % Normal case: Multiple slices, possibly multiple groups per slice
    for iSlice = 1:nSlices
        sliceIdx = outputSlices(iSlice);
        fprintf('\nProcessing slice %d/%d (data slice %d)...\n', iSlice, nSlices, sliceIdx);
        
        % Get coordinate information for this slice
        if length(projectionMatrixCoordinates_ARA) >= sliceIdx && ~isempty(projectionMatrixCoordinates_ARA{sliceIdx})
            coords = projectionMatrixCoordinates_ARA{sliceIdx};
            fprintf('  Found coordinates for slice %d\n', sliceIdx);
        else
            fprintf('  Warning: No coordinate information for slice %d\n', sliceIdx);
            continue;
        end
    end
end

% Process the coordinate information
if ~isempty(coords)
    fprintf('  Processing coordinates...\n');
        
        if length(coords) >= 3 && ~isempty(coords{3})
            ara_coordinate = coords{3}; % ARA level (this is already a slice number)
            atlas_slice_idx = round(ara_coordinate / 2); % Convert from 10um to 20um atlas
            fprintf('  ARA coordinate: %d (atlas slice at 20um: %d)\n', ara_coordinate, atlas_slice_idx);
            
            if atlas_slice_idx > 0 && atlas_slice_idx <= size(av_v2, 1)
                atlas_slice = squeeze(av_v2(atlas_slice_idx, :, :));
                fprintf('  Extracted atlas slice %d, size: %s\n', atlas_slice_idx, mat2str(size(atlas_slice)));
                
                % Get coordinate grids
                if length(coords) >= 2
                    x_coords = coords{1}; % ML coordinates
                    y_coords = coords{2}; % DV coordinates
                    fprintf('  X coords range: %.1f to %.1f (%d points)\n', min(x_coords), max(x_coords), length(x_coords));
                    fprintf('  Y coords range: %.1f to %.1f (%d points)\n', min(y_coords), max(y_coords), length(y_coords));
                    
                    % Process each group in this slice
                    for iGroup = 1:nGroups
                        fprintf('    Processing group %d/%d...\n', iGroup, nGroups);
                        
                        % Get slice data for this group
                        if nGroups > 1 && nActualSlices == 1
                            % Groups in 3rd dimension of 3D array
                            slice_data = projectionMatrix_array(:, :, iGroup);
                        elseif ndims(projectionMatrix_array) == 4
                            % 4D array: [pixels, pixels, groups, slices]
                            slice_data = projectionMatrix_array(:, :, iGroup, 1); % Use slice 1 for now
                        else
                            % 3D array with slices in 3rd dimension
                            slice_data = projectionMatrix_array(:, :, 1); % Use slice 1 for now
                        end
                        
                        fprintf('      Slice data size: %s, range: %.6f to %.6f\n', ...
                            mat2str(size(slice_data)), min(slice_data(:)), max(slice_data(:)));
                        
                        % Convert data coordinates from 10um to 20um space
                        x_coords_20um = x_coords / resolution_factor;
                        y_coords_20um = y_coords / resolution_factor;
                        
                        % Process each subregion for this group
                        fprintf('      Checking %d CP subregions in atlas slice...\n', nSubregions);
                        
                        % Check what subregion IDs are actually in this atlas slice
                        unique_ids_in_slice = unique(atlas_slice(:));
                        fprintf('      Atlas slice contains %d unique region IDs\n', length(unique_ids_in_slice));
                        
                        % Show some example IDs in the slice
                        fprintf('      First 10 IDs in slice: %s\n', mat2str(unique_ids_in_slice(1:min(10,end))));
                        
                        % Show some CP subregion IDs we're looking for
                        fprintf('      Looking for CP IDs: %s\n', mat2str(cp_subregions.id(1:min(5,end))')');
                        
                        nSubregionsWithVoxels = 0;
                        for iSubregion = 1:nSubregions
                            subregion_id = cp_subregions.id(iSubregion);
                            
                            % Find voxels belonging to this subregion in 20um atlas
                            subregion_mask = (atlas_slice == subregion_id);
                            
                            if any(subregion_mask(:))
                                nSubregionsWithVoxels = nSubregionsWithVoxels + 1;
                                if nSubregionsWithVoxels <= 5  % Debug first 5 subregions
                                    fprintf('      Found %d voxels for subregion %s (ID: %d)\n', ...
                                        sum(subregion_mask(:)), cp_subregions.acronym{iSubregion}, subregion_id);
                                end
                                
                                % Get subregion voxel coordinates in atlas space
                                [atlas_y, atlas_x] = find(subregion_mask);
                                
                                % Debug coordinate ranges
                                if nSubregionsWithVoxels == 1  % Only for first subregion
                                    fprintf('        Atlas voxel X range: %d to %d\n', min(atlas_x), max(atlas_x));
                                    fprintf('        Atlas voxel Y range: %d to %d\n', min(atlas_y), max(atlas_y));
                                    fprintf('        Data X coords (20um): %.1f to %.1f\n', min(x_coords_20um), max(x_coords_20um));
                                    fprintf('        Data Y coords (20um): %.1f to %.1f\n', min(y_coords_20um), max(y_coords_20um));
                                end
                                
                                % Initialize intensity collection
                                subregion_intensities = [];
                                
                                % For each atlas voxel in the subregion
                                nVoxelsInRange = 0;
                                nVoxelsWithData = 0;
                                for iVoxel = 1:length(atlas_x)
                                    atlas_x_phys = atlas_x(iVoxel);
                                    atlas_y_phys = atlas_y(iVoxel);
                                    
                                    % Find corresponding data grid location
                                    if atlas_x_phys >= min(x_coords_20um) && atlas_x_phys <= max(x_coords_20um) && ...
                                       atlas_y_phys >= min(y_coords_20um) && atlas_y_phys <= max(y_coords_20um)
                                        
                                        nVoxelsInRange = nVoxelsInRange + 1;
                                        
                                        % Interpolate to find data grid index
                                        data_x = interp1(x_coords_20um, 1:length(x_coords), atlas_x_phys, 'nearest', 'extrap');
                                        data_y = interp1(y_coords_20um, 1:length(y_coords), atlas_y_phys, 'nearest', 'extrap');
                                        
                                        % Ensure indices are within bounds
                                        data_x = round(max(1, min(size(slice_data, 1), data_x)));
                                        data_y = round(max(1, min(size(slice_data, 2), data_y)));
                                        
                                        % Extract intensity value
                                        intensity_val = slice_data(data_x, data_y);
                                        if ~isnan(intensity_val) && intensity_val > 0
                                            subregion_intensities = [subregion_intensities, intensity_val];
                                            nVoxelsWithData = nVoxelsWithData + 1;
                                        end
                                    end
                                end
                                
                                % Debug output for first subregion
                                if nSubregionsWithVoxels == 1
                                    fprintf('        Voxels in coordinate range: %d/%d\n', nVoxelsInRange, length(atlas_x));
                                    fprintf('        Voxels with data > 0: %d\n', nVoxelsWithData);
                                end
                                
                                if ~isempty(subregion_intensities)
                                    % Calculate mean intensity for this slice and group
                                    mean_intensity = mean(subregion_intensities);
                                    voxel_count = length(subregion_intensities);
                                    
                                    % Store results (using linear indexing for group combination)
                                    if nGroups > 1 && nActualSlices == 1
                                        % Groups in 3rd dimension: use group index directly
                                        result_idx = iGroup;
                                    else
                                        % Normal case: use group+slice combination
                                        result_idx = (iGroup-1)*nSlices + 1; % Use slice 1 for now
                                    end
                                    
                                    if result_idx <= size(subregionResults.mean_intensities, 2)
                                        subregionResults.mean_intensities(iSubregion, result_idx) = mean_intensity;
                                        subregionResults.voxel_counts(iSubregion, result_idx) = voxel_count;
                                    end
                                    
                                    % Accumulate for global calculations
                                    global_intensity_sums(iSubregion) = global_intensity_sums(iSubregion) + sum(subregion_intensities);
                                    global_voxel_counts(iSubregion) = global_voxel_counts(iSubregion) + voxel_count;
                                    
                                    fprintf('        %s (group %d): mean=%.4f, voxels=%d\n', ...
                                        cp_subregions.acronym{iSubregion}, iGroup, mean_intensity, voxel_count);
                                end
                            end
                        end
                    end
                else
                    fprintf('  Warning: Missing coordinate information\n');
                end
            else
                fprintf('  Warning: Atlas slice index %d out of bounds\n', atlas_slice_idx);
            end
        else
            fprintf('  Warning: Missing ARA coordinate\n');
        end
else
    fprintf('  Warning: No coordinate information available\n');
end

%% Calculate global averages
fprintf('\nCalculating global averages...\n');

for iSubregion = 1:nSubregions
    if global_voxel_counts(iSubregion) > 0
        global_mean = global_intensity_sums(iSubregion) / global_voxel_counts(iSubregion);
        globalResults.mean_intensities(iSubregion) = global_mean;
        globalResults.total_voxel_counts(iSubregion) = global_voxel_counts(iSubregion);
        
        fprintf('  %s (global): mean=%.4f, total_voxels=%d\n', ...
            cp_subregions.acronym{iSubregion}, global_mean, global_voxel_counts(iSubregion));
    end
end

%% Summary statistics
fprintf('\n=== SUMMARY ===\n');
fprintf('Slices analyzed: %d\n', nSlices);
fprintf('Groups analyzed: %d\n', nGroups);
fprintf('Subregions found: %d\n', nSubregions);

% Find subregions with data
subregions_with_data = sum(~isnan(subregionResults.mean_intensities), 2) > 0;
fprintf('Subregions with data: %d\n', sum(subregions_with_data));

if any(subregions_with_data)
    fprintf('\nTop 5 subregions by global mean intensity:\n');
    [~, sort_idx] = sort(globalResults.mean_intensities, 'descend', 'MissingPlacement', 'last');
    for i = 1:min(5, sum(~isnan(globalResults.mean_intensities)))
        idx = sort_idx(i);
        fprintf('  %d. %s: %.4f\n', i, cp_subregions.acronym{idx}, globalResults.mean_intensities(idx));
    end
end

%% Create bar plots for each group
fprintf('\nCreating bar plots for each group...\n');

% Find subregions with data for plotting
subregions_with_data_idx = find(sum(~isnan(subregionResults.mean_intensities), 2) > 0);
nSubregionsWithData = length(subregions_with_data_idx);

if nSubregionsWithData > 0
    % Create figure
    fig = figure('Name', 'CP Subregion Analysis by Group', 'Position', [100, 100, 1200, 800]);
    
    % Calculate subplot layout
    nCols = min(3, nGroups); % Max 3 columns
    nRows = ceil(nGroups / nCols);
    
    for iGroup = 1:nGroups
        subplot(nRows, nCols, iGroup);
        
        % Extract data for this group (assuming single slice for now)
        if nActualSlices == 1
            group_data = subregionResults.mean_intensities(subregions_with_data_idx, iGroup);
        else
            % For multiple slices, take mean across slices for this group
            group_cols = ((iGroup-1)*nSlices + 1):(iGroup*nSlices);
            group_data = nanmean(subregionResults.mean_intensities(subregions_with_data_idx, group_cols), 2);
        end
        
        % Get subregion names for plotting
        subregion_names = subregionResults.subregion_acronyms(subregions_with_data_idx);
        
        % Remove NaN values for plotting
        valid_idx = ~isnan(group_data);
        plot_data = group_data(valid_idx);
        plot_names = subregion_names(valid_idx);
        
        if ~isempty(plot_data)
            % Create bar plot
            bars = bar(plot_data, 'FaceColor', [0.2, 0.6, 0.8]);
            
            % Get group title - try to use region names if available
            if ~isempty(inputRegions) && ~isempty(regionGroups)
                % Find which input regions belong to this group
                if length(unique(regionGroups)) >= iGroup
                    uniqueGroups = unique(regionGroups);
                    currentGroupNum = uniqueGroups(iGroup);
                    regionsInGroup = find(regionGroups == currentGroupNum);
                    
                    if ~isempty(regionsInGroup)
                        regionNames = inputRegions(regionsInGroup);
                        if length(regionNames) == 1
                            group_title = regionNames{1};
                        else
                            group_title = strjoin(regionNames, '+');
                        end
                    else
                        group_title = sprintf('Group %d', iGroup);
                    end
                else
                    group_title = sprintf('Group %d', iGroup);
                end
            else
                group_title = sprintf('Group %d', iGroup);
            end
            
            % Customize plot
            title(group_title, 'FontWeight', 'bold', 'FontSize', 12);
            ylabel('Mean Fluorescence Intensity', 'FontSize', 10);
            xlabel('CP Subregions', 'FontSize', 10);
            
            % Set x-axis labels
            set(gca, 'XTick', 1:length(plot_names));
            set(gca, 'XTickLabel', plot_names);
            set(gca, 'XTickLabelRotation', 45);
            
            % Add value labels on bars
            for i = 1:length(plot_data)
                if plot_data(i) > 0
                    text(i, plot_data(i) + max(plot_data)*0.02, sprintf('%.3f', plot_data(i)), ...
                        'HorizontalAlignment', 'center', 'FontSize', 8);
                end
            end
            
            % Set consistent y-axis limits across all subplots
            if iGroup == 1
                all_data = subregionResults.mean_intensities(:);
                all_data = all_data(~isnan(all_data));
                if ~isempty(all_data)
                    global_max = max(all_data);
                    global_min = min(all_data);
                else
                    global_max = 1;
                    global_min = 0;
                end
                if global_max > global_min
                    ylim([0, global_max * 1.1]);
                end
            else
                if exist('global_max', 'var') && global_max > global_min
                    ylim([0, global_max * 1.1]);
                end
            end
            
            % Remove grid lines for cleaner look
        else
            % No data for this group
            text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 14);
            title(sprintf('Group %d', iGroup), 'FontWeight', 'bold', 'FontSize', 12);
            set(gca, 'XTick', [], 'YTick', []);
        end
    end
    
    % Add overall title
    sgtitle('Mean Fluorescence Intensity per CP Subregion by Group', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Adjust subplot spacing
    set(gcf, 'Color', 'white');
    
    fprintf('Bar plot created with %d groups and %d subregions with data\n', nGroups, nSubregionsWithData);
else
    fprintf('No subregions with data found - skipping plot creation\n');
end

fprintf('Analysis complete!\n\n');

end
