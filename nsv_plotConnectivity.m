function nsv_plotConnectivity(experimentData, allenAtlasPath, inputRegion, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

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
            xlabel('Left-to-right (a.u.)');
            hold on;
        elseif iChunk == round(numberOfChunks/2)
            ylabel('Posterior-to-anterior (a.u.)');
            hold on;
        end
    else
        subplot(1, numberOfChunks, iChunk)
        if iChunk == 1
            ylabel('Posterior-to-anterior (a.u.)');
            hold on;
        elseif iChunk == round(numberOfChunks/2)
            xlabel('Left-to-right (a.u.)');
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
    hold on;


    projection_view_lims(iChunk, 1, :) = xlim;
    projection_view_lims(iChunk, 2, :) = ylim;
    projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
        (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 1, 2), ...
        projection_view_lims(iChunk, 2, 1): ...
        (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 2, 2)};
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
    region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
        1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere
    % AP, DV, ML -> ML, AP, DV

    [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
        = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
    regionLocation(2, :) = regionLocation(2, :) + round(chunks_region(iChunk)) - 1;
    thisChunk_x = regionLocation(projection_views(1, 1), :);
    thisChunk_DV = regionLocation(projection_views(1, 2), :);

    subplot(1, numberOfChunks, iChunk)
    boundary_projection{iChunk} = boundary(thisChunk_x', ...
        thisChunk_DV', 0);

    plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
        regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
        'Color', plot_structure_color);
    axis equal
    axis square
    axis image

    projection_view_lims(iChunk, 1, :) = xlim;
    projection_view_lims(iChunk, 2, :) = ylim;
    projection_view_bins{iChunk} = {projection_view_lims(iChunk, 1, 1): ...
        (projection_view_lims(iChunk, 1, 2) - projection_view_lims(iChunk, 1, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 1, 2), ...
        projection_view_lims(iChunk, 2, 1): ...
        (projection_view_lims(iChunk, 2, 2) - projection_view_lims(iChunk, 2, 1)) / numberOfPixels: ...
        projection_view_lims(iChunk, 2, 2)};

end


theseLocations_ap_dv_ml = experimentData(:, :, 1:projectionGridSize(3)/2) + ...
    experimentData(:, :, projectionGridSize(3):-1:projectionGridSize(3)/2+1); % collapse ML so we only have one hemisphere

projectionMatrix = {};

% Loop over each chunk
for iChunk = 1:numberOfChunks
    nBinsX = numberOfPixels;
    nBinsY = numberOfPixels;
    
    projectionMatrix{iChunk} = zeros(nBinsX, nBinsY);

    % Retrieve bin edges for the current chunk
    xEdges = projection_view_bins{iChunk}{1};
    yEdges = projection_view_bins{iChunk}{2};

    thisdiff = mean(diff(chunks_region(:)));

    % find voxels that fit inside
    if ~any(xEdges == 0) && ~any(yEdges == 0)
        projtemp = permute(squeeze(nanmean(theseLocations_ap_dv_ml(round((chunks_region(iChunk) - thisdiff)./10):round((chunks_region(iChunk) + thisdiff)./10), ...
            round(yEdges/10), round(xEdges/10)), 1)), [2, 1]);

        projectionMatrix{iChunk} = projtemp;
    end

end

% normalize each slice (TO DO QQ)
for iChunk = 1:numberOfChunks
    % Extract the current chunk
    currentChunk = projectionMatrix{iChunk}(:, :);

    % Place the normalized chunk back into the matrix
    projectionMatrix{iChunk}(:, :) = currentChunk;
end

%% get and plot all fluorescence
figProjection = figure('Name', 'Fluorescence intensity');

for iChunk = 1:numberOfChunks
    thisCmap_limits = [-max(max(projectionMatrix{iChunk})), max(max(projectionMatrix{iChunk}))];

    binnedArrayPixelSmooth = projectionMatrix{iChunk};
    % remove any data points outside of the ROI
    clearvars regionLocation
    % get structure boundaries and plot outline
    region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
        1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]); % / 2 to only get one hemispehere
    % AP, DV, ML -> ML, AP, DV

    [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] ...
        = ind2sub(size(region_area), find(region_area)); %ML, AP, DV
    isIN = nan(size(binnedArrayPixelSmooth, 1), size(binnedArrayPixelSmooth, 2));
    for iPixelX = 1:size(binnedArrayPixelSmooth, 1)
        for iPixelY = 1:size(binnedArrayPixelSmooth, 2)
            isIN(iPixelX, iPixelY) = inpolygon(projection_view_bins{iChunk}{1}(iPixelX), ...
                projection_view_bins{iChunk}{2}(iPixelY), ...
                regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
                regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}));
        end
    end

    figure(figProjection);
    subplot(1, numberOfChunks, iChunk)
    binnedArrayPixelSmooth(isIN == 0) = NaN;
    ax = gca;
    ax.YColor = 'w'; % Red
    ax.XColor = 'w'; % Red
    im = imagesc(projection_view_bins{iChunk}{1}, projection_view_bins{iChunk}{2}, ...
        binnedArrayPixelSmooth');
    set(im, 'AlphaData', ~isnan(get(im, 'CData')));

    set(gca, 'color', [0.5, 0.5, 0.5]);

    originalColormap = brewermap([], '*RdBu');

    colormap(originalColormap);


    try
        caxis(thisCmap_limits)
    catch
        caxis([-1, 1])
    end
    hold on;

    %colorbar
    clearvars binnedArrayPixel
    hold on;
    boundary_projection{iChunk} = boundary(regionLocation(projection_views(iChunk, 1), :)', ...
        regionLocation(projection_views(iChunk, 2), :)', 0);
    plot(regionLocation(projection_views(iChunk, 1), boundary_projection{iChunk}), ...
        regionLocation(projection_views(iChunk, 2), boundary_projection{iChunk}), ...
        'Color', plot_structure_color, 'LineWidth', 2);

    axis equal
    axis square
    axis image
    ax.XLabel.Color = [0, 0, 0];
    ax.YLabel.Color = [0, 0, 0];
    nColors = numel(ax.YTickLabel);
    cm = [0, 0, 0];
    for i = 1:nColors
        ax.YTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.YTickLabel{i})];
    end

    nColors = numel(ax.XTickLabel);
    cm = [0, 0, 0];
    for i = 1:nColors
        ax.XTickLabel{i} = ['\color[rgb]', sprintf('{%f,%f,%f}%s', cm, ax.XTickLabel{i})];
    end

    clearvars isIN
    try
        caxis(thisCmap_limits)
    catch
        caxis([-1, 1])
    end
    set(gca, 'color', [0.5, 0.5, 0.5]);
    xlim([projection_view_bins{iChunk}{1}(1), projection_view_bins{iChunk}{1}(end)])
    ylim([projection_view_bins{iChunk}{2}(1), projection_view_bins{iChunk}{2}(end)])

    axis off;

    % Alternatively, to remove only specific components:
    set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
    xlabel(''); % Remove x-axis label
    ylabel(''); % Remove y-axis label
    thisAp_slice = round(nanmean(chunks_region(iChunk:iChunk+1)), 2);
    %  coordinates = ya_convert_allen_to_paxinos(values_toConvert, allenOrBrainreg)
    if iChunk == 1

        title(['Bregma: ', num2str(thisAp_slice)]); % Remove title
    elseif iChunk == numberOfChunks
        title([num2str(thisAp_slice), 'mm']); % Remove title
    else
        title([num2str(thisAp_slice)]); % Remove title
    end

end


% set x and ylims
xlims_region = nan(numberOfChunks, 2);
ylims_region = nan(numberOfChunks, 2);
for iChunk = 1:numberOfChunks
    subplot(size(1, 2), numberOfChunks, iChunk)
    xlims_region(iChunk, :) = xlim;
    ylims_region(iChunk, :) = ylim;
end

diff_xlims_region = diff(xlims_region');
diff_ylims_region = diff(ylims_region');
for iChunk = 1:numberOfChunks
    subplot(1, numberOfChunks, iChunk)
    xlims_here = (max(diff_xlims_region) - diff_xlims_region(iChunk)) ./ 2;
    xlim([xlims_region(iChunk, 1) - xlims_here, xlims_region(iChunk, 2) + xlims_here])

    ylims_here = (max(diff_ylims_region) - diff_ylims_region(iChunk)) ./ 2;
    ylim([ylims_region(iChunk, 1) - ylims_here, ylims_region(iChunk, 2) + ylims_here])
    if iChunk == 1
        axis_length_mm = 1;
        one_pixel_x = (diff(projection_view_bins{iChunk}{1}(2:3)));
        one_pixel_x_um = one_pixel_x * 25;
        one_pixel_y = (diff(projection_view_bins{iChunk}{2}(2:3)));
        one_pixel_y_um = one_pixel_y * 25;
        axis_length_atlas_units_x = (axis_length_mm * 1000) / (one_pixel_x_um);
        axis_length_atlas_units_y = (axis_length_mm * 1000) / (one_pixel_y_um);
        prettify_addScaleBars(axis_length_atlas_units_x, axis_length_atlas_units_y, ...
            [num2str(axis_length_mm), 'mm'], [num2str(axis_length_mm), 'mm'], 'topLeft', '', '')
    end

end


end
