function plotMultiRegionInjections(experimentImgs, allenAtlasPath, inputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color, normalizationMethod, experimentRegionInfo)
% plotMultiRegionInjections - Plot injection sites for multiple regions in one figure
%
% Creates a single figure with one row per input region showing injection sites

% Create main figure
mainFig = figure('Name', 'Multi-Region Injection Sites', 'Color', 'w');

nRegions = length(inputRegions);

fprintf('Creating combined injection plot for %d regions...\n', nRegions);

% Loop through each region and create subplot
for iRegion = 1:nRegions
    fprintf('Processing region %d/%d: %s\n', iRegion, nRegions, inputRegions{iRegion});
    
    % For each region, we need to call plotConnectivity but capture the data
    % instead of letting it create its own figure
    
    % Temporarily redirect figure creation
    oldFig = get(0, 'CurrentFigure');
    
    % Create temporary invisible figure
    tempFig = figure('Visible', 'off');
    
    try
        % Call plotConnectivity for this specific region
        bsv.plotConnectivity(experimentImgs, allenAtlasPath, inputRegions(iRegion), numberOfSlices, numberOfPixels, plane, ...
            regionOnly, smoothing, colorLimits, color, normalizationMethod, [], [], experimentRegionInfo, false);
        
        % Get all axes from the temporary figure
        tempAxes = findobj(tempFig, 'Type', 'axes');
        
        % Copy each subplot (chunk) to the main figure
        for iChunk = 1:numberOfSlices
            if length(tempAxes) >= iChunk
                % Calculate position in main figure: nRegions rows, numberOfSlices columns
                subplotPos = (iRegion - 1) * numberOfSlices + iChunk;
                
                % Switch to main figure and create subplot
                figure(mainFig);
                ax = subplot(nRegions, numberOfSlices, subplotPos);
                
                % Get the corresponding temporary axis (they're in reverse order from findobj)
                tempAx = tempAxes(end - iChunk + 1);
                
                % Copy all children (images, plots, etc.) from temp axis to main axis
                copyobj(get(tempAx, 'Children'), ax);
                
                % Copy axis properties
                set(ax, 'XLim', get(tempAx, 'XLim'), 'YLim', get(tempAx, 'YLim'));
                set(ax, 'CLim', get(tempAx, 'CLim'));
                colormap(ax, get(tempAx, 'Colormap'));
                
                % Set axis properties for clean display
                axis(ax, 'equal', 'tight', 'off');
                
                % Add labels
                if iChunk == 1
                    % First column: add region name
                    ylabel(sprintf('%s', inputRegions{iRegion}), 'FontWeight', 'bold', 'FontSize', 10);
                end
                
                if iRegion == 1 && ~isempty(get(tempAx, 'Title'))
                    % First row: copy title
                    title(get(get(tempAx, 'Title'), 'String'), 'FontSize', 9);
                end
            end
        end
        
    catch ME
        fprintf('Warning: Error processing region %s: %s\n', inputRegions{iRegion}, ME.message);
    end
    
    % Close temporary figure
    close(tempFig);
    
    % Restore original figure
    if ~isempty(oldFig)
        figure(oldFig);
    end
end

% Add overall title and colorbar
figure(mainFig);
sgtitle('Injection Sites by Region', 'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar to the last subplot
if nRegions > 0 && numberOfSlices > 0
    subplot(nRegions, numberOfSlices, nRegions * numberOfSlices);
    cb = colorbar('eastoutside');
    
    % Label colorbar based on normalization
    switch lower(normalizationMethod)
        case 'injectionvolume'
            labelText = 'Injection intensity (norm. by volume)';
        case 'injectionintensity'
            labelText = 'Injection intensity (norm. by intensity)';
        case 'none'
            labelText = 'Raw injection intensity';
        otherwise
            labelText = 'Injection intensity';
    end
    
    cb.Label.String = labelText;
    cb.Label.FontSize = 9;
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = 'bottom';
end

fprintf('Multi-region injection plot completed!\n');

end