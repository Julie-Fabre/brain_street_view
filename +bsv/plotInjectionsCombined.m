function plotInjectionsCombined(experimentImgs, allenAtlasPath, inputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationMethod, experimentRegionInfo, atlasType, atlasResolution)
% plotInjectionsCombined - Plot injection sites for multiple regions in ONE figure
%
% Creates a single figure with one row per input region showing injection sites
% This function builds the plot manually without calling plotConnectivity

% Handle optional atlas parameters
if nargin < 13 || isempty(atlasType)
    atlasType = 'allen'; % Default to Allen atlas
end
if nargin < 14 || isempty(atlasResolution)
    atlasResolution = 10; % Default to 10um resolution
end

%% Construct atlas filenames based on type and resolution
switch lower(atlasType)
    case 'allen'
        if atlasResolution == 10
            annotationFile = 'annotation_volume_10um_by_index.npy';
            structureFile = 'structure_tree_safe_2017.csv';
        elseif atlasResolution == 20
            annotationFile = 'annotation_volume_v2_20um_by_index.npy';
            structureFile = 'UnifiedAtlas_Label_ontology_v2.csv';
        else
            error('Unsupported Allen atlas resolution: %d. Supported resolutions are 10 and 20 um.', atlasResolution);
        end
    otherwise
        % For custom atlases, construct filename from type and resolution
        annotationFile = sprintf('%s_annotation_%dum.npy', atlasType, atlasResolution);
        structureFile = sprintf('%s_structure_tree.csv', atlasType);
end

%% Load atlas
av = readNPY([allenAtlasPath, filesep, annotationFile]);
st = loadStructureTree([allenAtlasPath, filesep, structureFile]);
atlas_slice_spacing = atlasResolution; % Use the provided resolution

%% Create single combined figure
mainFig = figure('Name', 'Combined Injection Sites - All Regions', 'Color', 'w');
nRegions = length(inputRegions);

fprintf('Creating ONE combined figure for %d regions...\n', nRegions);

%% Process each region and create subplots in the main figure
for iRegion = 1:nRegions
    fprintf('Adding region %d/%d: %s to combined plot\n', iRegion, nRegions, inputRegions{iRegion});
    
    % Find structure indices for this region
    curr_plot_structure_idx = find(contains(st.acronym, inputRegions{iRegion}));
    keepStruct = false(size(curr_plot_structure_idx, 1), 1);
    for iCurrStruct = 1:size(curr_plot_structure_idx, 1)
        curr_struct = st.acronym{curr_plot_structure_idx(iCurrStruct)};
        keepStruct(iCurrStruct) = strcmp(curr_struct(1:length(inputRegions{iRegion})), inputRegions{iRegion});
    end
    curr_plot_structure_idx = curr_plot_structure_idx(keepStruct);
    
    if isempty(curr_plot_structure_idx)
        fprintf('Warning: No structure found for region %s\n', inputRegions{iRegion});
        continue;
    end
    
    plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx(1)}, 2, [])') ./ 255;
    
    % Get chunk limits for this region
    structureLimits = find(ismember(av(:, :, 1:1:end/2), curr_plot_structure_idx));
    [APvalues, ~, MLvalues] = ind2sub(size(av), structureLimits);
    
    if strcmp(plane, 'coronal')
        curr_limits = [min(APvalues), max(APvalues)];
        chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfSlices:curr_limits(2);
    else
        curr_limits = [min(MLvalues), max(MLvalues)];
        chunks_region = curr_limits(1):(curr_limits(2) - curr_limits(1)) / numberOfSlices:curr_limits(2);
    end
    
    % Collapse hemispheres for data
    halfSlices = floor(size(experimentImgs, 3)/2);
    theseLocations_ap_dv_ml = experimentImgs(:, :, 1:halfSlices) + experimentImgs(:, :, end:-1:end-halfSlices+1);
    
    % Process each chunk (slice) for this region
    for iChunk = 1:numberOfSlices
        % Calculate subplot position: nRegions rows, numberOfSlices columns
        subplotPos = (iRegion - 1) * numberOfSlices + iChunk;
        subplot(nRegions, numberOfSlices, subplotPos);
        hold on;
        
        % Get structure boundaries for this chunk
        if strcmp(plane, 'coronal')
            region_area = permute(ismember(av(round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1)), ...
                1:1:end, 1:1:end/2), curr_plot_structure_idx), [3, 1, 2]);
        else
            region_area = permute(ismember(av(1:1:end, ...
                1:1:end, round(chunks_region(iChunk)):1:round(chunks_region(iChunk+1))), ...
                curr_plot_structure_idx), [3, 1, 2]);
        end
        
        [regionLocation(1, :), regionLocation(2, :), regionLocation(3, :)] = ind2sub(size(region_area), find(region_area));
        
        if strcmp(plane, 'coronal')
            regionLocation(2, :) = regionLocation(2, :) + round(chunks_region(iChunk)) - 1;
        else
            regionLocation(1, :) = regionLocation(1, :) + round(chunks_region(iChunk)) - 1;
        end
        
        % Get projection data for this chunk
        thisdiff = mean(diff(chunks_region(:)));
        
        if strcmp(plane, 'coronal')
            % Create bins for this chunk
            x_range = [min(regionLocation(1, :)), max(regionLocation(1, :))];
            y_range = [min(regionLocation(3, :)), max(regionLocation(3, :))];
            
            if diff(x_range) > 0 && diff(y_range) > 0
                xEdges = linspace(x_range(1), x_range(2), numberOfPixels+1);
                yEdges = linspace(y_range(1), y_range(2), numberOfPixels+1);
                
                % Extract and average data
                dataSlice = theseLocations_ap_dv_ml(round((chunks_region(iChunk) - thisdiff)/10):round((chunks_region(iChunk) + thisdiff)/10), ...
                    round(yEdges/10), round(xEdges/10));
                meanData = nanmean(dataSlice, 1);
                projtemp = permute(squeeze(meanData), [2, 1]);
                
                % Create mask for region boundary
                [X, Y] = meshgrid(xEdges(1:end-1), yEdges(1:end-1));
                if size(regionLocation, 2) > 2
                    boundary_idx = boundary(regionLocation(1, :)', regionLocation(3, :)', 0);
                    mask = inpolygon(X, Y, regionLocation(1, boundary_idx), regionLocation(3, boundary_idx));
                    projtemp(~mask') = NaN;
                end
                
                % Plot the data
                imagesc(xEdges(1:end-1), yEdges(1:end-1), projtemp');
                
                % Plot region boundary
                if size(regionLocation, 2) > 2
                    plot(regionLocation(1, boundary_idx), regionLocation(3, boundary_idx), ...
                        'Color', plot_structure_color, 'LineWidth', 2);
                end
                
                % Set colormap and limits
                colormap(gca, flipud(gray(256)));
                set(gca, 'color', [0.5, 0.5, 0.5]);
                clim([0, 1]);
                
                % Make NaN values grey
                set(gca, 'color', [0.5, 0.5, 0.5]);
            end
        end
        
        % Format axes
        axis equal; axis tight; axis off;
        if strcmp(plane, 'coronal')
            set(gca, 'YDir', 'reverse');
        end
        
        % Add labels
        if iChunk == 1
            ylabel(sprintf('%s', inputRegions{iRegion}), 'FontWeight', 'bold', 'FontSize', 10, 'Rotation', 90);
        end
        
        if iRegion == 1
            slice_ARA = round(mean(chunks_region(iChunk:iChunk+1))/10);
            title(sprintf('ARA %d', slice_ARA), 'FontSize', 9);
        end
    end
end

% Add overall title and colorbar
sgtitle('Injection Sites by Region', 'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar
subplot(nRegions, numberOfSlices, nRegions * numberOfSlices);
cb = colorbar('eastoutside');
cb.Label.String = 'Injection Intensity';
cb.Label.FontSize = 9;
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = 'bottom';

fprintf('Combined injection plot completed in ONE figure!\n');

end