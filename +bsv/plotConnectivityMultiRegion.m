function [projectionMatrix_array, projectionMatrixCoordinates_ARA] = plotConnectivityMultiRegion(experimentData, allenAtlasPath, outputRegions, numberOfChunks, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationInfo, inputRegions, regionGroups, experimentRegionInfo, normalizeByGroup, customSlices, sliceAveraging)
% plotConnectivityMultiRegion - Plot multiple brain regions in one figure, one row per region
%
% This function creates a combined plot with each output region in its own row.
% It's a wrapper around plotConnectivity that handles multiple regions.

if nargin < 17 || isempty(sliceAveraging)
    sliceAveraging = 0;
end
if nargin < 16 || isempty(customSlices)
    customSlices = [];
end
if nargin < 15 || isempty(normalizeByGroup)
    normalizeByGroup = false;
end

% Ensure outputRegions is a cell array
if ischar(outputRegions) || isstring(outputRegions)
    outputRegions = {char(outputRegions)};
end

nOutputRegions = length(outputRegions);

% Create a new figure for the combined plot
figProjection = figure('Name', 'Multi-Region Fluorescence Intensity', 'Color', 'w');

% Initialize output arrays
projectionMatrix_array = {};
projectionMatrixCoordinates_ARA = {};

fprintf('Plotting %d output regions in combined view...\n', nOutputRegions);

% Loop through each output region
for iOutputRegion = 1:nOutputRegions
    fprintf('Processing region %d/%d: %s\n', iOutputRegion, nOutputRegions, outputRegions{iOutputRegion});
    
    % Get data for this specific region by calling plotConnectivity
    % We'll capture the figure and extract the relevant plots
    
    % Create temporary figure for this region
    tempFig = figure('Visible', 'off');
    
    try
        [tempProjMatrix, tempProjCoords] = bsv.plotConnectivity(experimentData, allenAtlasPath, ...
            outputRegions(iOutputRegion), numberOfChunks, numberOfPixels, plane, regionOnly, ...
            smoothing, colorLimits, color, normalizationInfo, inputRegions, regionGroups, ...
            experimentRegionInfo, normalizeByGroup, customSlices, sliceAveraging);
        
        % Store the results
        projectionMatrix_array{iOutputRegion} = tempProjMatrix;
        projectionMatrixCoordinates_ARA{iOutputRegion} = tempProjCoords;
        
        % Get the data from the temporary figure
        tempAxes = findall(tempFig, 'Type', 'axes');
        
        % Copy each subplot to the main figure
        for iChunk = 1:numberOfChunks
            if ~isempty(tempAxes) && length(tempAxes) >= iChunk
                % Calculate subplot position in combined figure
                % nOutputRegions rows, numberOfChunks columns
                subplotIdx = (iOutputRegion - 1) * numberOfChunks + iChunk;
                
                % Create subplot in main figure
                figure(figProjection);
                subplot(nOutputRegions, numberOfChunks, subplotIdx);
                
                % Copy content from temporary figure
                if length(tempAxes) >= iChunk
                    tempAx = tempAxes(end - iChunk + 1); % Reverse order due to findall
                    
                    % Copy images and plot elements
                    copyobj(allchild(tempAx), gca);
                    
                    % Copy axis properties
                    xlim(tempAx.XLim);
                    ylim(tempAx.YLim);
                    set(gca, 'XTick', tempAx.XTick, 'YTick', tempAx.YTick);
                    set(gca, 'XTickLabel', tempAx.XTickLabel, 'YTickLabel', tempAx.YTickLabel);
                    colormap(gca, tempAx.Colormap);
                    clim(tempAx.CLim);
                    
                    % Set title and labels
                    if iOutputRegion == 1
                        title(tempAx.Title.String);
                    end
                    
                    if iChunk == 1
                        ylabel(sprintf('%s', outputRegions{iOutputRegion}), 'FontWeight', 'bold');
                    end
                    
                    % Remove ticks and labels for cleaner look
                    set(gca, 'XTick', [], 'YTick', []);
                    axis off;
                end
            end
        end
        
    catch ME
        fprintf('Error processing region %s: %s\n', outputRegions{iOutputRegion}, ME.message);
    end
    
    % Close temporary figure
    close(tempFig);
end

% Add colorbar to the combined figure
if nOutputRegions > 0 && numberOfChunks > 0
    % Add colorbar to the last subplot
    subplot(nOutputRegions, numberOfChunks, nOutputRegions * numberOfChunks);
    cb = colorbar('eastoutside');
    
    % Create informative label based on normalization method
    if normalizeByGroup
        labelText = ['Projection intensity' newline '(normalized by group)'];
    else
        switch lower(normalizationInfo)
            case 'injectionvolume'
                labelText = ['Projection intensity' newline '(norm. by injection volume)'];
            case 'injectionintensity'
                labelText = ['Projection intensity' newline '(norm. by injection intensity)'];
            case 'none'
                labelText = 'Raw projection intensity';
            otherwise
                labelText = 'Projection intensity';
        end
    end
    
    cb.Label.String = labelText;
    cb.Label.FontSize = 9;
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = 'bottom';
end

fprintf('Multi-region plot completed with %d regions\n', nOutputRegions);

end