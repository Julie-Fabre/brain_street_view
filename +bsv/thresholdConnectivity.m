function [projectionMatrix_array, projectionMatrixCoordinates_ARA] = thresholdConnectivity(experimentData, allenAtlasPath, ...
    inputRegion, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, threshold, thresholdMethod, normalizationMethod, dataFetchNormalization)
% Enhanced thresholding and connectivity visualization
%
% INPUTS:
%   threshold - threshold value OR percentile (0-100) for percentile method
%   thresholdMethod - 'absolute' (default), 'percentile', 'zscore', 'relative'
%   normalizationMethod - 'none' (default), 'region', 'zscore', 'robust'
%
% THRESHOLDING METHODS:
%   'absolute' - Direct threshold value
%   'percentile' - Use percentile of data (threshold = percentile value 0-100)
%   'zscore' - Threshold at N standard deviations above mean (threshold = N)
%   'relative' - Threshold as percentage of maximum (threshold = 0-1)
%
% NORMALIZATION METHODS:
%   'none' - No additional normalization
%   'region' - Scale each region to [0,1] independently  
%   'zscore' - Z-score normalization (mean=0, std=1)
%   'robust' - Robust scaling using median and IQR

% Set defaults for new parameters
if nargin < 12 || isempty(thresholdMethod)
    thresholdMethod = 'absolute';
end

if nargin < 13 || isempty(normalizationMethod)
    normalizationMethod = 'none';
end

% Handle optional dataFetchNormalization parameter
if nargin < 14 || isempty(dataFetchNormalization)
    dataFetchNormalization = 'unknown';
end

%% load Allen atlas
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
atlas_slice_spacing = 10; % 10 um/ slice

%% current region info
curr_plot_structure_idx = find(contains(st.acronym, inputRegion));
% check
keepStruct = false(size(curr_plot_structure_idx, 1), 1);
for iCurrStruct = 1:size(curr_plot_structure_idx, 1)
    curr_struct = st.acronym{curr_plot_structure_idx(iCurrStruct)};
    keepStruct(iCurrStruct) = strcmp(curr_struct(1:size(inputRegion{:}, 2)), inputRegion);
end
curr_plot_structure_idx = curr_plot_structure_idx(keepStruct);

%curr_plot_structure = st.id(curr_plot_structure_idx);
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx(1)}, 2, [])') ./ 255;

%% Get chunk limits in AP (if coronal) or ML (if sagital)
structureLimits = find(ismember(av(:, :, 1:1:end/2), curr_plot_structure_idx));
[APvalues, ~, MLvalues] = ind2sub(size(av), structureLimits); % ap x dv x ml
if strcmp(plane, 'coronal')
    curr_limits = [min(APvalues), max(APvalues)];
    chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
    projection_views = repmat([1, 2], numberOfChunks, 1); % ML x AP
elseif strcmp(plane, 'sagital')
    curr_limits = [min(MLvalues), max(MLvalues)];
    chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
    projection_views = repmat([2, 1], numberOfChunks, 1); % AP x ML
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

if size(experimentData, 4) == 1
    theseLocations_ap_dv_ml = zeros(size(experimentData)./[1, 1, 2]);
else
    theseLocations_ap_dv_ml = zeros(size(experimentData)./[1, 1, 2, 1]);
end
for iGroup = 1:size(experimentData, 4)
    theseLocations_ap_dv_ml(:, :, :, iGroup) = ...
        experimentData(:, :, 1:projectionGridSize(3)/2, iGroup) + ...
        experimentData(:, :, projectionGridSize(3):-1:projectionGridSize(3)/2+1, iGroup); % collapse ML so we only have one hemisphere
end

projectionMatrix = {};

% Loop over each chunk
for iChunk = 1:numberOfChunks
    nBinsX = numberOfPixels+1;
    nBinsY = numberOfPixels+1;

    projectionMatrix{iChunk} = zeros(nBinsX, nBinsY);

    % Retrieve bin edges for the current chunk
    xEdges = projection_view_bins{iChunk}{1};
    yEdges = projection_view_bins{iChunk}{2};

    thisdiff = mean(diff(chunks_region(:)));

    % find voxels that fit inside
    if ~any(xEdges == 0) && ~any(yEdges == 0)
        if strcmp(plane, 'coronal')
            projtemp = permute(squeeze(nanmean( ...
                theseLocations_ap_dv_ml(round( ...
                (chunks_region(iChunk) - thisdiff)./10):round((chunks_region(iChunk) + thisdiff)./10), ...
                round(yEdges/10),...
                round(xEdges/10), :), 1)), [2, 1, 3]); % AP x DV x ML ->  DV x AP x ML
        else
            projtemp = permute(squeeze(nanmean( ...
                theseLocations_ap_dv_ml(round(xEdges/10),...
                round(yEdges/10), ...
                round((chunks_region(iChunk) - thisdiff)./10):...
                round((chunks_region(iChunk) + thisdiff)./10), :), 3)), [2, 1, 3]);

        end
        projectionMatrix{iChunk} = projtemp;
    end

end

%% Apply enhanced normalization and thresholding
fprintf('\n🔧 PROCESSING PARAMETERS\n');
fprintf('Threshold method: %s\n', thresholdMethod);
fprintf('Normalization method: %s\n', normalizationMethod);

% Apply normalization first
for iChunk = 1:numberOfChunks
    for iGroup = 1:size(experimentData, 4)
        currentChunk = projectionMatrix{iChunk}(:, :, iGroup);
        
        switch lower(normalizationMethod)
            case 'region'
                % Scale to [0,1] independently
                minVal = min(currentChunk(:));
                maxVal = max(currentChunk(:));
                if maxVal > minVal
                    currentChunk = (currentChunk - minVal) / (maxVal - minVal);
                end
                
            case 'zscore'
                % Z-score normalization
                meanVal = mean(currentChunk(:), 'omitnan');
                stdVal = std(currentChunk(:), 'omitnan');
                if stdVal > 0
                    currentChunk = (currentChunk - meanVal) / stdVal;
                end
                
            case 'robust'
                % Robust scaling using median and IQR
                medianVal = median(currentChunk(:), 'omitnan');
                q75 = prctile(currentChunk(:), 75);
                q25 = prctile(currentChunk(:), 25);
                iqr = q75 - q25;
                if iqr > 0
                    currentChunk = (currentChunk - medianVal) / iqr;
                end
                
            case 'none'
                % No additional normalization
                
            otherwise
                warning('Unknown normalization method: %s. Using none.', normalizationMethod);
        end
        
        projectionMatrix{iChunk}(:, :, iGroup) = currentChunk;
    end
end

% Calculate adaptive threshold based on method
adaptiveThreshold = threshold; % default
allData = [];

% Collect all data for threshold calculation
for iChunk = 1:numberOfChunks
    for iGroup = 1:size(experimentData, 4)
        chunkData = projectionMatrix{iChunk}(:, :, iGroup);
        allData = [allData; chunkData(:)];
    end
end

% Remove NaN values for threshold calculation
allData = allData(~isnan(allData));

switch lower(thresholdMethod)
    case 'percentile'
        adaptiveThreshold = prctile(allData, threshold);
        fprintf('Calculated %gth percentile threshold: %.4f\n', threshold, adaptiveThreshold);
        
    case 'zscore'
        meanVal = mean(allData);
        stdVal = std(allData);
        adaptiveThreshold = meanVal + threshold * stdVal;
        fprintf('Calculated %g-sigma threshold: %.4f (mean=%.4f, std=%.4f)\n', ...
            threshold, adaptiveThreshold, meanVal, stdVal);
        
    case 'relative'
        maxVal = max(allData);
        adaptiveThreshold = threshold * maxVal;
        fprintf('Calculated relative threshold: %.4f (%g%% of max=%.4f)\n', ...
            adaptiveThreshold, threshold*100, maxVal);
        
    case 'absolute'
        fprintf('Using absolute threshold: %.4f\n', threshold);
        
    otherwise
        warning('Unknown threshold method: %s. Using absolute.', thresholdMethod);
end

% Calculate and display processing statistics
numThresholdedVoxels = sum(allData > adaptiveThreshold);
totalVoxels = length(allData);
percentThresholded = (numThresholdedVoxels / totalVoxels) * 100;

fprintf('Data statistics:\n');
fprintf('  • Total voxels: %d\n', totalVoxels);
fprintf('  • Voxels above threshold: %d (%.1f%%)\n', numThresholdedVoxels, percentThresholded);
fprintf('  • Data range: [%.4f, %.4f]\n', min(allData), max(allData));
fprintf('  • Mean ± SD: %.4f ± %.4f\n', mean(allData), std(allData));
fprintf('─────────────────────────────────\n');

%% get and plot all fluorescence
figProjection = figure('Name', 'Thresholded Connectivity', 'Color', 'w');
figOriginal = figure('Name', 'Original with Threshold Boundary', 'Color', 'w');

nGroups = size(experimentData, 4);

% Store original data for the second figure
originalProjectionMatrix = projectionMatrix;

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

    % Create threshold mask for polygon overlay using adaptive threshold
    thresholdMask = originalProjectionMatrix{iChunk}(:, :, 1) > adaptiveThreshold;
    
    % Apply adaptive threshold
    projectionMatrix{iChunk}(:, :, iGroup) = projectionMatrix{iChunk}(:, :, iGroup) > adaptiveThreshold;

    for iGroup = 1:nGroups
        % setup plotting axis
        figure(figProjection);
        subplot(nGroups, numberOfChunks, (iGroup - 1)*numberOfChunks+iChunk)
        hold on;
        ax = gca;

        % Enhanced colormap and scaling for thresholded data
        thisCmap_limits = [0, 1]; % Binary for thresholded data

        % remove any data points outside of the region
        binnedArrayPixel = projectionMatrix{iChunk}(:, :, iGroup);
        binnedArrayPixel(isIN == 0) = NaN;
        [isIN_i, isIN_j] = find(isIN==0);
        
        % Set corresponding elements in projectionMatrix to 0
        projectionMatrix{iChunk}(sub2ind(size(projectionMatrix{iChunk}(:,:,iGroup)), isIN_i, isIN_j, repmat(iGroup, size(isIN_i)))) = 0;

        % smooth - QQ TODO
        binnedArrayPixelSmooth = binnedArrayPixel;

        % plot data
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            binnedArrayPixelSmooth');

        % make NaN values (outside of ROI) grey
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);

        % Enhanced colormap for thresholded data - black/white only
        colormap(gca, flipud(gray(256))); % Black = high, white = low (flipped)
        clim(thisCmap_limits);
        
        % Add colorbar with proper labels for thresholded figure
        if iChunk == numberOfChunks && iGroup == nGroups
            cb = colorbar('eastoutside');
            
            % Create informative label based on normalization (same as plotConnectivity)
            switch lower(dataFetchNormalization)
                case 'injectionvolume'
                    labelText = ['Projection intensity' newline '(norm. by injection' newline 'volume)'];
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
            
            cb.Label.String = labelText;
            cb.Label.FontSize = 9;
            cb.Label.Rotation = 270;
            cb.Label.VerticalAlignment = 'bottom';
            
            % Position colorbar to not overlap with plots
            cb.Position(1) = cb.Position(1) + 0.02;
        end

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
        xlabel(''); % Remove y-axis label
        set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
        axis off;

        % title
        this_slice_ARA = round(nanmean(chunks_region(iChunk:iChunk+1))./10); %  coordinates = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
        if iChunk == 1
            if strcmp(plane, 'coronal')
                title(['ARA level (cor.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
            else
                title(['ARA level (sag.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
            end
        elseif iChunk == numberOfChunks
            title([num2str(this_slice_ARA)]);
        else
            title([num2str(this_slice_ARA)]);
        end

        sliceARAs(iChunk) = this_slice_ARA;

        clearvars binnedArrayPixelSmooth binnedArrayPixel

    end
    
    %% Plot original data with threshold boundary overlay
    for iGroup = 1:nGroups
        % setup plotting axis for original figure
        figure(figOriginal);
        subplot(nGroups, numberOfChunks, (iGroup - 1)*numberOfChunks+iChunk)
        hold on;
        ax = gca;

        % Enhanced colormap limits for original data - consistent across regions
        if strcmp(normalizationMethod, 'none')
            % Use global max across all data for consistency
            globalMax = max(allData);
            thisCmap_limits = [0, globalMax];
        else
            % Normalized data - use consistent scale
            if strcmp(normalizationMethod, 'region')
                thisCmap_limits = [0, 1];
            elseif strcmp(normalizationMethod, 'zscore')
                thisCmap_limits = [-3, 3]; % ±3 sigma range
            else % robust
                thisCmap_limits = [-3, 3]; % ±3 IQR range
            end
        end

        % remove any data points outside of the region
        binnedArrayPixel_orig = originalProjectionMatrix{iChunk}(:, :, iGroup);
        binnedArrayPixel_orig(isIN == 0) = NaN;

        % smooth - QQ TODO
        binnedArrayPixelSmooth_orig = binnedArrayPixel_orig;

        % plot original data
        im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
            binnedArrayPixelSmooth_orig');

        % make NaN values (outside of ROI) grey
        set(im, 'AlphaData', ~isnan(get(im, 'CData')));
        set(gca, 'color', [0.5, 0.5, 0.5]);

        % Enhanced colormap for original data - black/white only
        colormap(gca, flipud(gray(256))); % Black = high, white = low (flipped)
        clim(thisCmap_limits);
        
        % Add colorbar with proper labels for original figure
        if iChunk == numberOfChunks && iGroup == nGroups
            cb = colorbar('eastoutside');
            
            % Create informative label based on normalization (same as plotConnectivity)
            switch lower(dataFetchNormalization)
                case 'injectionvolume'
                    labelText = ['Projection intensity' newline '(norm. by injection' newline 'volume)'];
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
            
            cb.Label.String = labelText;
            cb.Label.FontSize = 9;
            cb.Label.Rotation = 270;
            cb.Label.VerticalAlignment = 'bottom';
            
            % Position colorbar to not overlap with plots
            cb.Position(1) = cb.Position(1) + 0.02;
        end

        % plot region boundaries
        plot(regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
            regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}), ...
            'Color', plot_structure_color, 'LineWidth', 2);

        % Plot threshold boundary as ONE dotted red polygon following square edges
        if any(thresholdMask(:))
            % Get bin coordinates and widths
            xBins = projection_view_bins{iChunk}{1};
            yBins = projection_view_bins{iChunk}{2};
            
            if length(xBins) > 1 && length(yBins) > 1
                xBinWidth = xBins(2) - xBins(1);
                yBinWidth = yBins(2) - yBins(1);
                
                % Create a slightly larger grid to handle edges properly
                expandedMask = padarray(thresholdMask, [1 1], 0, 'both');
                
                % Use bwboundaries to trace the boundary on the pixel grid
                boundaries = bwboundaries(expandedMask, 'noholes');
                
                if ~isempty(boundaries)
                    % Get the largest boundary (main region)
                    boundarySizes = cellfun(@length, boundaries);
                    [~, maxIdx] = max(boundarySizes);
                    regionBoundary = boundaries{maxIdx};
                    
                    % Convert boundary indices back to coordinate space
                    % regionBoundary contains [row, col] indices in the expanded mask
                    % Convert to actual x,y coordinates
                    xCoords = [];
                    yCoords = [];
                    
                    for i = 1:size(regionBoundary, 1)
                        % Get expanded mask indices (1-based)
                        expandedRow = regionBoundary(i, 1);
                        expandedCol = regionBoundary(i, 2);
                        
                        % Convert back to original mask indices (accounting for padding)
                        origRow = expandedRow - 1;
                        origCol = expandedCol - 1;
                        
                        % Convert to coordinate space with pixel edges
                        if origRow >= 1 && origRow <= length(xBins) && ...
                           origCol >= 1 && origCol <= length(yBins)
                            % Get pixel center
                            xCenter = xBins(origRow);
                            yCenter = yBins(origCol);
                            
                            % Calculate boundary coordinate (on pixel edge)
                            xCoord = xCenter + (expandedCol - origCol - 1.5) * xBinWidth/2;
                            yCoord = yCenter + (expandedRow - origRow - 1.5) * yBinWidth/2;
                        else
                            % Handle edge cases
                            if origRow < 1
                                xCoord = xBins(1) - xBinWidth/2;
                            elseif origRow > length(xBins)
                                xCoord = xBins(end) + xBinWidth/2;
                            else
                                xCoord = xBins(origRow);
                            end
                            
                            if origCol < 1
                                yCoord = yBins(1) - yBinWidth/2;
                            elseif origCol > length(yBins)
                                yCoord = yBins(end) + yBinWidth/2;
                            else
                                yCoord = yBins(origCol);
                            end
                        end
                        
                        xCoords = [xCoords; xCoord];
                        yCoords = [yCoords; yCoord];
                    end
                    
                    % Plot the single boundary polygon
                    if length(xCoords) > 2
                        plot([xCoords; xCoords(1)], [yCoords; yCoords(1)], ...
                             'r:', 'LineWidth', 2);
                    end
                end
            end
        end

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
        xlabel(''); % Remove y-axis label
        set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
        axis off;

        % title
        this_slice_ARA = round(nanmean(chunks_region(iChunk:iChunk+1))./10);
        if iChunk == 1
            if strcmp(plane, 'coronal')
                title(['Original - ARA level (cor.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
            else
                title(['Original - ARA level (sag.): ', num2str(this_slice_ARA)]); % Allen Reference Atlas level
            end
        elseif iChunk == numberOfChunks
            title([num2str(this_slice_ARA)]);
        else
            title([num2str(this_slice_ARA)]);
        end

        clearvars binnedArrayPixelSmooth_orig binnedArrayPixel_orig

    end

end


% set x and ylims for both figures
xlims_region = nan(numberOfChunks, 2);
ylims_region = nan(numberOfChunks, 2);

% Get limits from thresholded figure
figure(figProjection);
for iChunk = 1:numberOfChunks
    subplot(nGroups, numberOfChunks, (iGroup - 1)*numberOfChunks+iChunk)
    xlims_region(iChunk, :) = xlim;
    ylims_region(iChunk, :) = ylim;
end

diff_xlims_region = diff(xlims_region');
diff_ylims_region = diff(ylims_region');

% Apply limits to both figures
for figHandle = [figProjection, figOriginal]
    figure(figHandle);
    for iChunk = 1:numberOfChunks
        for iGroup = 1:nGroups
            subplot(nGroups, numberOfChunks, (iGroup - 1)*numberOfChunks+iChunk)

            xlims_here = (max(diff_xlims_region) - diff_xlims_region(iChunk)) ./ 2;
            xlim([xlims_region(iChunk, 1) - xlims_here, xlims_region(iChunk, 2) + xlims_here])

            ylims_here = (max(diff_ylims_region) - diff_ylims_region(iChunk)) ./ 2;
            ylim([ylims_region(iChunk, 1) - ylims_here, ylims_region(iChunk, 2) + ylims_here])

            shrink_factor_x = diff(xlims_region(iChunk, :)) ./ (diff(xlims_region(iChunk, :)) + xlims_here * 2);
            %shrink_factor_y = diff(ylims_region(iChunk, :)) ./ (diff(ylims_region(iChunk, :)) + ylims_here*2);
            if iChunk == 1 && iGroup == 1 && figHandle == figProjection
                axis_length_mm = 1;
                one_pixel_x = (diff(projection_view_bins{iChunk}{1}(2:3)));
                one_pixel_x_um = one_pixel_x / 2.5 / shrink_factor_x; % QQ why 2.5 here
                % one_pixel_y = (diff(projection_view_bins{iChunk}{2}(2:3)));
                % one_pixel_y_um = one_pixel_y;
                axis_length_atlas_units_x = (axis_length_mm * 1000) / (one_pixel_x_um);
                % axis_length_atlas_units_y = (axis_length_mm * 1000) / (one_pixel_y_um);
                prettify_addScaleBars(axis_length_atlas_units_x, 0, ...
                    [num2str(axis_length_mm), 'mm'], '', 'topLeft', '', '')
            end
        end
    end
end
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
