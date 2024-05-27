function nsv_plotConnectivity(experimentData, allenAtlasPath, inputRegion, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color)

%% load Allen atlas
[tv, av, st, bregma] = ya_loadAllenAtlas(allenAtlasPath); % QQ bregma wrong in brainreg atlas

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

curr_plot_structure = st.id(curr_plot_structure_idx);
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx(1)}, 2, [])') ./ 255;

%% Get chunk limits

structureLimits = find(ismember(av(:, :, 1:1:end/2), curr_plot_structure_idx));
[APvalues, DVvalues, MLvalues] = ind2sub(size(av), structureLimits); % ap x dv x ml
if strcmp(plane, 'coronal')
    curr_limits = [min(APvalues), max(APvalues)];
    chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
    projection_views = repmat([1, 2], numberOfChunks, 1); % AP x ML
elseif strcmp(plane, 'sagital')
    curr_limits = [min(MLvalues), max(MLvalues)];
    chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfChunks:curr_limits(2);
    projection_views = repmat([2, 1], numberOfChunks, 1); % AP x ML
end

% initialize variables
boundary_projection = cell(3, 1);
projection_view_bins = cell(3, 1);
projection_view_lims = nan(3, 2, 2);

figOutline1 = figure('Name', 'Chunk AP/ML limits');
xlims_region = nan(numberOfChunks, 2);
ylims_region = nan(numberOfChunks, 2);

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
            xlabel('Left-to-right (a.u.)'); hold on;
        elseif iChunk == round(numberOfChunks/2)
            ylabel('Posterior-to-anterior (a.u.)'); 
            hold on;
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
        else
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            hold on;
        end
    else
        subplot(1, numberOfChunks, iChunk)
         if iChunk == 1
            ylabel('Posterior-to-anterior (a.u.)'); hold on;
        elseif iChunk == round(numberOfChunks/2)
            xlabel('Left-to-right (a.u.)'); 
            hold on;
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
        else
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            hold on;
        end
    end
    boundary_projection{iChunk} = boundary(thisChunk_x', ...
        thisChunk_y', 0);
 if strcmp(plane, 'coronal')
    plot(regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
        -regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
        'Color', plot_structure_color);
 else
     plot(regionLocation(projection_views(1, 2), boundary_projection{iChunk}), ...
        regionLocation(projection_views(1, 1), boundary_projection{iChunk}), ...
        'Color', plot_structure_color);
 end

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


figOutline2 = figure();

figProjection = figure();

end
