function [projectionMatrix_array, projectionMatrixCoordinates_ARA] = plotConnectivity(experimentData, allenAtlasPath, outputRegion,...
    numberOfChunks, numberOfPixels, plane,...
    regionOnly, smoothing, colorLimits, color, normalizationInfo,...
    inputRegions, regionGroups, experimentRegionInfo, normalizeByGroup, customSlices, sliceAveraging)
% plotConnectivity - Plot brain connectivity data with optional region grouping
%
% Parameters:
%   outputRegion - Single target region to visualize (e.g., 'CP' for striatum)
%   inputRegions - Cell array of source region names (e.g., {'SS', 'MOp', 'MOs'}) 
%   regionGroups (optional) - Numeric array indicating which group each input region belongs to.
%                           Length must match the number of input regions.
%                           For example, with 4 regions: [1, 2, 2, 3] means:
%                           - Region 1 projections in row 1
%                           - Regions 2&3 projections averaged in row 2 
%                           - Region 4 projections in row 3
%                           If empty or not provided, all input regions are combined into one group.
%   experimentRegionInfo - Structure containing experiment region information
%   normalizeByGroup - Boolean to normalize within groups (default: false)
%   customSlices (optional) - Array of specific Allen atlas slice numbers to plot (e.g., [40, 45, 50])
%                           If empty or not provided, uses automatic chunking based on numberOfChunks
%   sliceAveraging (optional) - Number of slices to average around each specified slice (default: 0)
%                             For example, 2 means average from slice-2 to slice+2
%
% this function needs cleaning up + commenting

% Handle optional normalization info parameter
if nargin < 11 || isempty(normalizationInfo)
    normalizationInfo = 'unknown';
end

% Handle optional input regions parameter  
if nargin < 12 || isempty(inputRegions)
    inputRegions = []; % No input regions specified
end

% Handle optional region groups parameter
if nargin < 13 || isempty(regionGroups)
    regionGroups = []; % No grouping, use original behavior
end

% Handle optional experiment region info parameter
if nargin < 14 || isempty(experimentRegionInfo)
    experimentRegionInfo = []; % No experiment region info provided
end

% Handle optional normalizeByGroup parameter
if nargin < 15 || isempty(normalizeByGroup)
    normalizeByGroup = false; % Default: don't normalize by group
end

% Handle optional customSlices parameter
if nargin < 16 || isempty(customSlices)
    customSlices = []; % Use automatic chunking
end

% Handle optional sliceAveraging parameter  
if nargin < 17 || isempty(sliceAveraging)
    sliceAveraging = 0; % No averaging by default
end

%% load Allen atlas
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
atlas_slice_spacing = 10; % 10 um/ slice

%% Handle output region (target region for visualization)
if ischar(outputRegion) || isstring(outputRegion)
    outputRegion = {char(outputRegion)};
end

% Find structure indices for output region (defines the plotting area)
curr_plot_structure_idx = find(contains(st.acronym, outputRegion));
% check for exact match
keepStruct = false(size(curr_plot_structure_idx, 1), 1);
for iCurrStruct = 1:size(curr_plot_structure_idx, 1)
    curr_struct = st.acronym{curr_plot_structure_idx(iCurrStruct)};
    keepStruct(iCurrStruct) = strcmp(curr_struct(1:length(outputRegion{1})), outputRegion{1});
end
curr_plot_structure_idx = curr_plot_structure_idx(keepStruct);
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx(1)}, 2, [])') ./ 255;

%% Handle input regions (source regions) and grouping
if ~isempty(inputRegions)
    % Handle input regions - ensure it's a cell array
    if ischar(inputRegions) || isstring(inputRegions)
        inputRegions = {char(inputRegions)};
    end

    % Handle region grouping
    nRegions = length(inputRegions);
    
    if ~isempty(regionGroups)
        % regionGroups is a numeric array indicating which group each region belongs to
        if length(regionGroups) ~= nRegions
            error('Length of regionGroups (%d) must match number of input regions (%d)', length(regionGroups), nRegions);
        end
        
        % Validate group numbers
        if min(regionGroups) < 1 || any(mod(regionGroups, 1) ~= 0)
            error('Region group numbers must be positive integers');
        end
        
        % Convert to cell array format for internal processing
        uniqueGroups = unique(regionGroups);
        nRegionGroups = length(uniqueGroups);
        regionGroupsCellArray = cell(nRegionGroups, 1);
        
        for i = 1:nRegionGroups
            regionGroupsCellArray{i} = find(regionGroups == uniqueGroups(i));
        end
        regionGroups = regionGroupsCellArray;
    else
        % No grouping, each region gets its own row (original behavior)
        nRegionGroups = nRegions;
        regionGroups = num2cell(1:nRegions);
    end
    
    % Debug output
    fprintf('Number of input regions: %d\n', nRegions);
    fprintf('Number of region groups: %d\n', nRegionGroups);
    for i = 1:nRegionGroups
        fprintf('Group %d contains regions: %s\n', i, mat2str(regionGroups{i}));
    end
else
    % No input regions specified - use original behavior (single row)
    nRegionGroups = 1;
    regionGroups = {1};
    fprintf('No input regions specified - using original behavior\n');
end

%% Get chunk limits in AP (if coronal) or ML (if sagital)
structureLimits = find(ismember(av(:, :, 1:1:end/2), curr_plot_structure_idx));
[APvalues, ~, MLvalues] = ind2sub(size(av), structureLimits); % ap x dv x ml

% Check if custom slices are provided
if ~isempty(customSlices)
    % Use custom slices instead of automatic chunking
    % Convert Allen atlas slice numbers to array indices (multiply by 10 for 10um spacing)
    customSlices_indices = customSlices * atlas_slice_spacing;
    
    % Create chunks_region with boundaries for each custom slice
    % Each slice will be centered with sliceAveraging on each side
    chunks_region = [];
    for i = 1:length(customSlices_indices)
        slice_start = customSlices_indices(i) - sliceAveraging * atlas_slice_spacing;
        slice_end = customSlices_indices(i) + sliceAveraging * atlas_slice_spacing;
        chunks_region = [chunks_region, slice_start, slice_end];
    end
    
    % Remove duplicates and sort
    chunks_region = unique(chunks_region);
    
    % Update numberOfChunks to match custom slices
    numberOfChunks = length(customSlices);
    
    % Set projection views based on plane
    if strcmp(plane, 'coronal')
        projection_views = repmat([1, 2], numberOfChunks, 1); % ML x AP
    elseif strcmp(plane, 'sagital')
        projection_views = repmat([2, 1], numberOfChunks, 1); % AP x ML
    end
    
    fprintf('Using custom slices: %s (with averaging of %d slices on each side)\n', mat2str(customSlices), sliceAveraging);
else
    % Original automatic chunking behavior
    if strcmp(plane, 'coronal')
        curr_limits = [min(APvalues), max(APvalues)];
        chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
        projection_views = repmat([1, 2], numberOfChunks, 1); % ML x AP
    elseif strcmp(plane, 'sagital')
        curr_limits = [min(MLvalues), max(MLvalues)];
        chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
        projection_views = repmat([2, 1], numberOfChunks, 1); % AP x ML
    end
end

% initialize variables
boundary_projection = cell(3, 1);
projection_view_bins = cell(3, 1);
projection_view_lims = nan(3, 2, 2);

if strcmp(plane, 'coronal')
    figOutline1 = figure('Name', 'Chunk AP limits');
else
    figOutline1 = figure('Name', 'Chunk ML limits');
end

for iChunk = 1:numberOfChunks
    clearvars regionLocation
    % get structure boundaries and plot outline
    if strcmp(plane, 'coronal')
        region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
            1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    else
        region_area = permute(ismember(av(1:1:end, ...
            1:1:end, round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1))), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    end
    % AP, DV, ML -> ML, AP, DV

    [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
        = ind2sub(size(region_area), find(region_area)); %ML, AP, DV

    thisChunk_x = regionLocation(projection_views(1, 1), :);
    thisChunk_y = regionLocation(projection_views(1, 2), :);

    if strcmp(plane, 'coronal')
        subplot(numberOfChunks, 1, iChunk)
        if iChunk == numberOfChunks
            xlabel('Lateral-to-medial (a.u.)');
            hold on;
        elseif iChunk == round(numberOfChunks/2)
            ylabel('Posterior-to-anterior (a.u.)');
            hold on;
        end
    else
        subplot(1, numberOfChunks, iChunk)
        if iChunk == 1
            ylabel('Anterior-to-posterior (a.u.)');
            hold on;
        elseif iChunk == round(numberOfChunks/2)
            xlabel('Lateral-to-medial (a.u.)');
            hold on;
        end
    end
    boundary_projection{iChunk} = boundary(thisChunk_x', ...
        thisChunk_y', 0);
    if strcmp(plane, 'coronal')
        plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
            regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
            'Color', plot_structure_color);
    else
        plot(regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
            regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
            'Color', plot_structure_color);
    end

    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    projection_view_lims(iChunk, 1, :) = xlim;
    projection_view_lims(iChunk, 2, :) = ylim;
    projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
        (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 1, 2), ...
        projection_view_lims(iChunk, 2, 1): ...
        (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 2, 2)};

    if strcmp(plane, 'coronal')
        set(gca, 'YDir', 'reverse');
    end
end

if strcmp(plane, 'coronal')
    prettify_plot('XLimits', 'all');
else
    prettify_plot('YLimits', 'all');
end

%% for each chunk, get ML x DV values (if coronal) or AP x DV values (if sagital)
projectionGridSize = [132, 80, 114]; % experiment data (from Allen) is in AP (100 um) * DV (100 um) * ML(100 um)
% debugging :
% figure();
% iSlice=1
% iSlice = iSlice +1
% imagesc(squeeze(experimentData(iSlice,:,:)))
if strcmp(plane, 'coronal')
    figOutline2 = figure('Name', 'Chunk ML x DV limits');
    projection_views = repmat([1, 3], numberOfChunks, 1); % ML x DV
else
    figOutline2 = figure('Name', 'Chunk AP x DV limits');
    projection_views = repmat([2, 3], numberOfChunks, 1); % AP x DV
end

% initialize variables
boundary_projection = cell(3, 1);
projection_view_bins = cell(3, 1);
projection_view_lims = nan(3, 2, 2);

for iChunk = 1:numberOfChunks
    clearvars regionLocation
    % get structure boundaries and plot outline
    if strcmp(plane, 'coronal')
        region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
            1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    else
        region_area = permute(ismember(av(1:1:end, ...
            1:1:end, round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1))), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    end

    [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
        = ind2sub(size(region_area), find(region_area)); % ML, AP, DV
    if strcmp(plane, 'coronal')
        regionLocation(2, :) = regionLocation(2, :) + round(chunks_region(iChunk)) - 1;
    else
        regionLocation(1, :) = regionLocation(1, :) + round(chunks_region(iChunk)) - 1;
    end


    thisChunk_x = regionLocation(projection_views(1, 1), :);
    thisChunk_DV = regionLocation(projection_views(1, 2), :);

    subplot(1, numberOfChunks, iChunk)
    boundary_projection{iChunk} = boundary(thisChunk_x', ...
        thisChunk_DV', 0);
    hold on;

    if iChunk == 1
        ylabel('Ventral-to-dorsal (a.u.)');
        if strcmp(plane, 'coronal')
            xlabel('Lateral-to-medial (a.u.)');
        else
            xlabel('Anterior-to-posterior (a.u.)');
        end
    elseif iChunk == round(numberOfChunks/2)
        if strcmp(plane, 'coronal')
            title('<- Anterior to posterior slices ->');
        else
            title('<- Lateral to medial slices ->');
        end          
    end


    plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
        regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
        'Color', plot_structure_color);

    projection_view_lims(iChunk, 1, :) = xlim;
    projection_view_lims(iChunk, 2, :) = ylim;
    projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
        (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 1, 2), ...
        projection_view_lims(iChunk, 2, 1): ...
        (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 2, 2)};

    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    hold on;

end

prettify_plot;

% Now experimentData should be [AP, DV, ML, regionGroups] if grouping was applied

if size(experimentData, 4) > 1
    % 4D data with region groups - collapse ML hemispheres for each group
    nGroups = size(experimentData, 4);
    theseLocations_ap_dv_ml = zeros(size(experimentData, 1), size(experimentData, 2), floor(size(experimentData, 3)/2), nGroups);
    
    halfSlices = floor(size(experimentData, 3)/2);
    
    for iGroup = 1:nGroups
        leftHemi = experimentData(:, :, 1:halfSlices, iGroup);
        rightHemi = experimentData(:, :, end:-1:end-halfSlices+1, iGroup);
        
        
        theseLocations_ap_dv_ml(:, :, :, iGroup) = leftHemi + rightHemi;
    end
else
    % Original 3D data - collapse hemispheres
    halfSlices = floor(size(experimentData, 3)/2);
    theseLocations_ap_dv_ml = zeros(size(experimentData, 1), size(experimentData, 2), halfSlices);
    theseLocations_ap_dv_ml = ...
        experimentData(:, :, 1:halfSlices) + ...
        experimentData(:, :, end:-1:end-halfSlices+1);
end


projectionMatrix = {};

% Loop over each chunk
for iChunk = 1:numberOfChunks
    nBinsX = numberOfPixels+1;
    nBinsY = numberOfPixels+1;

    % Initialize with proper dimensions to accommodate different source regions
    if size(theseLocations_ap_dv_ml, 4) > 1
        projectionMatrix{iChunk} = zeros(nBinsX, nBinsY, size(theseLocations_ap_dv_ml, 4));
    else
        projectionMatrix{iChunk} = zeros(nBinsX, nBinsY);
    end

    % Retrieve bin edges for the current chunk
    xEdges = projection_view_bins{iChunk}{1};
    yEdges = projection_view_bins{iChunk}{2};

    thisdiff = mean(diff(chunks_region(:)));

    % find voxels that fit inside
    if ~any(xEdges == 0) && ~any(yEdges == 0)
        if iChunk == 1
        end
        
        if strcmp(plane, 'coronal')
            % Extract data slice
            dataSlice = theseLocations_ap_dv_ml(round( ...
                (chunks_region(iChunk) - thisdiff)./10):round((chunks_region(iChunk) + thisdiff)./10), ...
                round(yEdges/10),...
                round(xEdges/10), :);
            if iChunk == 1
            end
            
            % Average over AP dimension but keep the 4th dimension (different regions)
            meanData = nanmean(dataSlice, 1);
            if iChunk == 1
            end
            
            % Squeeze and permute, but be careful to preserve the region dimension
            projtemp = permute(squeeze(meanData), [2, 1, 3]); % DV x ML x regions ->  ML x DV x regions
            if iChunk == 1
            end
        else
            % Similar for sagittal
            dataSlice = theseLocations_ap_dv_ml(round(xEdges/10),...
                round(yEdges/10), ...
                round((chunks_region(iChunk) - thisdiff)./10):...
                round((chunks_region(iChunk) + thisdiff)./10), :);
            meanData = nanmean(dataSlice, 3);
            projtemp = permute(squeeze(meanData), [2, 1, 3]);
        end
        projectionMatrix{iChunk} = projtemp;
    end

end

% normalize each slice (TO DO QQ)
for iChunk = 1:numberOfChunks
    % Keep the 3rd dimension for different source regions
    % Currently no normalization applied - preserve original structure
    % projectionMatrix{iChunk} already contains the correct data from projtemp
end

% Apply group-wise normalization if requested
if normalizeByGroup && size(projectionMatrix{1}, 3) > 1
    fprintf('Applying group-wise normalization...\n');
    
    % Find max value for each group across all chunks
    nGroups = size(projectionMatrix{1}, 3);
    groupMaxValues = zeros(nGroups, 1);
    
    for iGroup = 1:nGroups
        maxVal = 0;
        for iChunk = 1:numberOfChunks
            groupData = projectionMatrix{iChunk}(:, :, iGroup);
            maxVal = max(maxVal, max(groupData(:)));
        end
        groupMaxValues(iGroup) = maxVal;
    end
    
    % Normalize each group to 0-1 range
    for iChunk = 1:numberOfChunks
        for iGroup = 1:nGroups
            if groupMaxValues(iGroup) > 0
                projectionMatrix{iChunk}(:, :, iGroup) = ...
                    projectionMatrix{iChunk}(:, :, iGroup) / groupMaxValues(iGroup);
            end
        end
    end
end

%% get and plot all fluorescence
figProjection = figure('Name', 'Fluorescence intensity', 'Color', 'w');

% Get nGroups from experimentData for compatibility with return value logic
nGroups = size(experimentData, 4);
for iChunk = 1:numberOfChunks

    clearvars regionLocation isIN

    % get structure boundaries
    if strcmp(plane, 'coronal')
        region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
            1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    else
        region_area = permute(ismember(av(1:1:end, ...
            1:1:end, round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1))), ...
            curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere, ML x AP x DV
    end
    [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
        = ind2sub(size(region_area), find(region_area)); %ML, AP, DV


    isIN = nan(size(projectionMatrix{iChunk}, 1), size(projectionMatrix{iChunk}, 2));
    for iPixelX = 1:size(projectionMatrix{iChunk}, 1)
        for iPixelY = 1:size(projectionMatrix{iChunk}, 2)
            isIN(iPixelX, iPixelY) = inpolygon(projection_view_bins{iChunk}{1}(iPixelX), ...
                projection_view_bins{iChunk}{2}(iPixelY), ...
                regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
                regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}));
        end
    end

    for iRegionGroup = 1:nRegionGroups
        % Get the regions that belong to this region group
        regionsInThisGroup = regionGroups{iRegionGroup};
        
        if iChunk == 1  % Only print debug for first chunk to avoid spam
        end
        
        % Extract data for this region group (already pre-grouped in fetchConnectivityData)
        if size(projectionMatrix{iChunk}, 3) >= iRegionGroup
            avgData = projectionMatrix{iChunk}(:, :, iRegionGroup);
        else
            % Fallback to 2D data
            avgData = projectionMatrix{iChunk}(:, :);
        end
        
        % setup plotting axis
        figure(figProjection);
        subplot(nRegionGroups, numberOfChunks, (iRegionGroup - 1)*numberOfChunks+iChunk)
        hold on;
        ax = gca;

        % colormap limits
        maxValue = max(avgData(:));
        thisCmap_limits = [0, 1];

        % remove any data points outside of the region
        binnedArrayPixel = avgData;
        binnedArrayPixel(isIN == 0) = NaN;

        % smooth - QQ TODO
        binnedArrayPixelSmooth = binnedArrayPixel;

        % plot data
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            binnedArrayPixelSmooth');

        % make NaN values (outside of ROI) grey
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);

        % set colormap - always black and white (flipped: black = high, white = low)
        colormap(gca, flipud(gray(256)));

        % set colormap limits
        clim(thisCmap_limits)

        % plot region boundaries
        plot(regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
            regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}), ...
            'Color', plot_structure_color, 'LineWidth', 2);

        % prettify axis
        axis equal
        axis square
        axis image
        ax.XLabel.Color = [0, 0, 0];
        ax.YLabel.Color = [0, 0, 0];

        % yaxis
        if strcmp(plane, 'coronal')
            set(ax, 'YDir', 'reverse');
        end
        nColors = numel(ax.YTickLabel);
        for i = 1:nColors
            ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', [0, 0, 0], ax.YTickLabel{i})];
        end
        ylim([projection_view_bins{iChunk}{2}(1), projection_view_bins{iChunk}{2}(end)])
        ylabel(''); % Remove y-axis label

        % xaxis
        nColors = numel(ax.XTickLabel);
        for i = 1:nColors
            ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', [0, 0, 0], ax.XTickLabel{i})];
        end
        xlim([projection_view_bins{iChunk}{1}(1), projection_view_bins{iChunk}{1}(end)])
        xlabel(''); % Remove x-axis label
        set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
        axis off;

        % title - only show ARA level at the top row, group labels on the left
        if ~isempty(customSlices)
            % For custom slices, use the provided slice number directly
            this_slice_ARA = customSlices(iChunk);
        else
            % Original calculation for automatic chunking
            this_slice_ARA = round(nanmean(chunks_region(iChunk:iChunk+1))./10); %  coordinates = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
        end
        
        if iRegionGroup == 1
            % Top row: show ARA levels
            if iChunk == 1
                if strcmp(plane, 'coronal')
                    title(['ARA level (cor.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
                else
                    title(['ARA level (sag.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
                end
            else
                title([num2str(this_slice_ARA)]);
            end
        end
        
        sliceARAs(iChunk) = this_slice_ARA;
        
        % Add group labels on the left side
        if iChunk == 1 && ~isempty(inputRegions)
            % Left column: show group labels
            regionsInThisGroup = regionGroups{iRegionGroup};
            groupRegionNames = inputRegions(regionsInThisGroup);
            if length(groupRegionNames) == 1
                groupLabel = groupRegionNames{1};
            else
                % Split into two lines if too many regions
                if length(groupRegionNames) > 3
                    mid = ceil(length(groupRegionNames) / 2);
                    line1 = strjoin(groupRegionNames(1:mid), '+');
                    line2 = strjoin(groupRegionNames(mid+1:end), '+');
                    groupLabel = [line1 newline line2];
                else
                    groupLabel = strjoin(groupRegionNames, '+');
                end
            end
            
            % Add text annotation to the left of the plot
            xlims = xlim;
            ylims = ylim;
            text(xlims(1) - 0.15*(xlims(2)-xlims(1)), mean(ylims), groupLabel, ...
                'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', 'Rotation', 90, 'Color', 'black');
        end
        
        % Add colorbar for the last subplot in each row
        if iChunk == numberOfChunks && iRegionGroup == nRegionGroups
            cb = colorbar('eastoutside');
            
            % Create informative label based on normalization method and data characteristics
            % Try to auto-detect from data dimensions
            dataDims = size(experimentData);
            if length(dataDims) >= 4 && dataDims(4) > 1
                nExperiments = dataDims(4);
            else
                nExperiments = 1; % Default for combined/averaged data
            end
            
            % Create experiment text
            if nExperiments == 1
                experimentText = '1 exp';
            else
                experimentText = sprintf('%d exps', nExperiments);
            end
            
            % Determine appropriate label based on normalization
            if normalizeByGroup
                % Group normalization overrides other normalization for display
                labelText = ['Projection intensity' newline '(normalized by group)'];
            else
                switch lower(normalizationInfo)
                    case 'injectionvolume'
                        labelText = ['Projection intensity (' newline '(norm. by injection' newline 'volume)'];
                    case 'injectionintensity'
                        labelText = ['Projection intensity' newline '(norm. by injection' newline 'intensity)'];
                    case 'none'
                        labelText = 'Raw projection intensity';
                    otherwise
                        % Auto-detect based on data range
                        dataMax = max(experimentData(:));
                        if dataMax <= 1.1 && dataMax > 0.1
                            labelText = ['Projection intensity' newline '(normalized)'];
                        else
                            labelText = 'Projection intensity';
                        end
                end
            end
            
            cb.Label.String = labelText;
            cb.Label.FontSize = 9;
            cb.Label.Rotation = 270;
            cb.Label.VerticalAlignment = 'bottom';
            
            % Position colorbar to not overlap with plots
            cb.Position(1) = cb.Position(1) + 0.02; % Move slightly right
        end

        clearvars binnedArrayPixelSmooth binnedArrayPixel

    end


end


% set x and ylims
xlims_region = nan(numberOfChunks, 2);
ylims_region = nan(numberOfChunks, 2);
for iChunk = 1:numberOfChunks
    subplot(nRegionGroups, numberOfChunks, (nRegionGroups - 1)*numberOfChunks+iChunk)
    xlims_region(iChunk, :) = xlim;
    ylims_region(iChunk, :) = ylim;
end

diff_xlims_region = diff(xlims_region');
diff_ylims_region = diff(ylims_region');
for iChunk = 1:numberOfChunks
    for iRegionGroup = 1:nRegionGroups
        subplot(nRegionGroups, numberOfChunks, (iRegionGroup - 1)*numberOfChunks+iChunk)

        xlims_here = (max(diff_xlims_region) - diff_xlims_region(iChunk)) ./ 2;
        xlim([xlims_region(iChunk, 1) - xlims_here, xlims_region(iChunk, 2) + xlims_here])

        ylims_here = (max(diff_ylims_region) - diff_ylims_region(iChunk)) ./ 2;
        ylim([ylims_region(iChunk, 1) - ylims_here, ylims_region(iChunk, 2) + ylims_here])

        shrink_factor_x = diff(xlims_region(iChunk, :)) ./ (diff(xlims_region(iChunk, :)) + xlims_here * 2);
        %shrink_factor_y = diff(ylims_region(iChunk, :)) ./ (diff(ylims_region(iChunk, :)) + ylims_here*2);
        if iChunk == 1 && iRegionGroup == 1
            axis_length_mm = 1;
            % Scale bar needs to be in data coordinates (atlas voxel units)
            % 1mm = 1000 μm, and each atlas voxel = 10 μm
            % So 1mm = 100 atlas voxels
            axis_length_atlas_units_x = (axis_length_mm * 1000) / atlas_slice_spacing;
            
            prettify_addScaleBars(axis_length_atlas_units_x, 0, ...
                [num2str(axis_length_mm), 'mm'], '', 'topLeft', '', '')
        end
    end

end
% Return matrix - note that when grouping is used, the 4th dimension will reflect 
% the averaged values across groups within each region group, not the original groups
if nGroups == 1
    projectionMatrix_array = zeros(size(projectionMatrix{1},1), size(projectionMatrix{1},2), size(projectionMatrix,2));
    for iSlice = 1:size(projectionMatrix,2)
        projectionMatrix_array(:,:,iSlice) = cell2mat(projectionMatrix(iSlice));
    end
else
    projectionMatrix_array = zeros(size(projectionMatrix{1},1), size(projectionMatrix{1},2), size(projectionMatrix{1},3), size(projectionMatrix,2));
    for iSlice = 1:size(projectionMatrix,2)
        projectionMatrix_array(:,:,:,iSlice) = cell2mat(projectionMatrix(iSlice));
    end
end
projectionMatrixCoordinates_ARA = projection_view_bins;
for iSlice = 1:size(projectionMatrix,2)
    projectionMatrixCoordinates_ARA{iSlice}{3} = sliceARAs(iSlice)*10;
end
end
