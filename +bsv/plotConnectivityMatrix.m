function [connectivityMatrix, figHandle] = plotConnectivityMatrix(projectionData, allenAtlasPath, targetRegions, varargin)
% plotConnectivityMatrix - Create an N x M connectivity matrix heatmap
%
% Build a matrix showing connectivity strength from N source regions to
% M target regions and display it as a heatmap.
%
% Parameters:
%   projectionData - Structure with fields named by source region acronyms.
%                    Each field contains the output from fetchConnectivityData:
%                    projectionData.VISam = {combinedProjection, combinedInfo, individual, expInfo}
%                    or just the projection array directly.
%   allenAtlasPath - Path to the directory containing atlas annotation volume
%                    and structure tree files.
%   targetRegions  - Cell array of target region acronyms (e.g., {'CP', 'ACB', 'SNr'}).
%
% Optional Parameters (name-value pairs):
%   'metric'         - Connectivity metric: 'mean' (default), 'max', or 'volume'.
%   'atlasType'      - Atlas type (default 'allen').
%   'atlasResolution'- Atlas resolution in micrometres (default 10).
%   'normalizeRows'  - If true, normalize each row to [0, 1] (default false).
%   'normalizeCols'  - If true, normalize each column to [0, 1] (default false).
%   'colormap'       - Colormap name or matrix (default 'viridis').
%   'annotate'       - If true, display numeric values in each cell (default true).
%   'title'          - Plot title (default empty).
%
% Returns:
%   connectivityMatrix - 2D array of shape (nSources, nTargets) with connectivity values.
%   figHandle          - Handle to the figure.
%
% Example:
%   % Fetch data for source regions
%   sourceRegions = {'VISam', 'VISp'};
%   projectionData = struct();
%   for i = 1:length(sourceRegions)
%       src = sourceRegions{i};
%       expIds = bsv.findConnectivityExperiments({src});
%       [proj, ~] = bsv.fetchConnectivityData(expIds, saveLocation, '', ...
%           'injectionIntensity', false, '', allenAtlasPath);
%       projectionData.(src) = proj;
%   end
%
%   % Create connectivity matrix
%   targetRegions = {'CP', 'ACB', 'SNr'};
%   [matrix, fig] = bsv.plotConnectivityMatrix(projectionData, allenAtlasPath, ...
%       targetRegions, 'metric', 'mean', 'annotate', true);

%% Parse inputs
p = inputParser;
addRequired(p, 'projectionData', @isstruct);
addRequired(p, 'allenAtlasPath', @ischar);
addRequired(p, 'targetRegions', @iscell);
addParameter(p, 'metric', 'mean', @(x) ismember(x, {'mean', 'max', 'volume'}));
addParameter(p, 'atlasType', 'allen', @ischar);
addParameter(p, 'atlasResolution', 10, @isnumeric);
addParameter(p, 'normalizeRows', false, @islogical);
addParameter(p, 'normalizeCols', false, @islogical);
addParameter(p, 'colormap', 'viridis', @(x) ischar(x) || isnumeric(x));
addParameter(p, 'annotate', true, @islogical);
addParameter(p, 'title', '', @ischar);
parse(p, projectionData, allenAtlasPath, targetRegions, varargin{:});

metric = p.Results.metric;
atlasType = p.Results.atlasType;
atlasResolution = p.Results.atlasResolution;
normalizeRows = p.Results.normalizeRows;
normalizeCols = p.Results.normalizeCols;
cmapName = p.Results.colormap;
annotate = p.Results.annotate;
plotTitle = p.Results.title;

%% Load atlas
switch lower(atlasType)
    case 'allen'
        if atlasResolution == 10
            annotationFile = 'annotation_volume_10um_by_index.npy';
            structureFile = 'structure_tree_safe_2017.csv';
        elseif atlasResolution == 20
            annotationFile = 'annotation_volume_v2_20um_by_index.npy';
            structureFile = 'UnifiedAtlas_Label_ontology_v2.csv';
        else
            error('Unsupported Allen atlas resolution: %d', atlasResolution);
        end
    otherwise
        annotationFile = sprintf('%s_annotation_%dum.npy', atlasType, atlasResolution);
        structureFile = sprintf('%s_structure_tree.csv', atlasType);
end

av = readNPY([allenAtlasPath, filesep, annotationFile]);
st = loadStructureTree([allenAtlasPath, filesep, structureFile]);

%% Downsample atlas to projection grid resolution (100 um)
projResolution = 100; % um
scaleFactor = atlasResolution / projResolution;

if scaleFactor ~= 1.0
    % Downsample using nearest-neighbor interpolation
    newSize = round(size(av) * scaleFactor);
    [X, Y, Z] = ndgrid(1:size(av,1), 1:size(av,2), 1:size(av,3));
    [Xq, Yq, Zq] = ndgrid(linspace(1, size(av,1), newSize(1)), ...
                          linspace(1, size(av,2), newSize(2)), ...
                          linspace(1, size(av,3), newSize(3)));
    avDownsampled = interp3(Y, X, Z, double(av), Yq, Xq, Zq, 'nearest');
else
    avDownsampled = av;
end

%% Build masks for each target region
nTargets = length(targetRegions);
targetMasks = cell(nTargets, 1);

for t = 1:nTargets
    region = targetRegions{t};
    % Find structure indices (prefix match)
    matchIdx = find(contains(st.acronym, region));
    keepIdx = [];
    for i = 1:length(matchIdx)
        acronym = st.acronym{matchIdx(i)};
        if length(acronym) >= length(region) && strcmp(acronym(1:length(region)), region)
            keepIdx = [keepIdx; matchIdx(i)]; %#ok<AGROW>
        end
    end

    if isempty(keepIdx)
        warning('No structures found matching ''%s''', region);
        targetMasks{t} = false(size(avDownsampled));
    else
        targetMasks{t} = ismember(avDownsampled, keepIdx);
    end
end

%% Get source regions and build connectivity matrix
sourceRegions = fieldnames(projectionData);
nSources = length(sourceRegions);

% Voxel volume in mm^3 (100 um = 0.1 mm per side)
voxelVolume = (projResolution / 1000)^3;

connectivityMatrix = zeros(nSources, nTargets);

for i = 1:nSources
    src = sourceRegions{i};
    data = projectionData.(src);

    % Handle both cell array (from fetchConnectivityData) and direct array
    if iscell(data)
        proj = data{1}; % combined_projection
    else
        proj = data;
    end

    % Average across groups if multiple groups exist
    if ndims(proj) == 4
        proj = mean(proj, 4, 'omitnan');
    end

    for j = 1:nTargets
        mask = targetMasks{j};

        % Handle shape mismatch between projection and mask
        if ~isequal(size(proj), size(mask))
            % Resize mask to match projection
            [Xm, Ym, Zm] = ndgrid(1:size(mask,1), 1:size(mask,2), 1:size(mask,3));
            [Xq, Yq, Zq] = ndgrid(linspace(1, size(mask,1), size(proj,1)), ...
                                  linspace(1, size(mask,2), size(proj,2)), ...
                                  linspace(1, size(mask,3), size(proj,3)));
            maskResized = interp3(Ym, Xm, Zm, double(mask), Yq, Xq, Zq, 'nearest') > 0.5;
        else
            maskResized = mask;
        end

        maskedValues = proj(maskResized);
        maskedValues = maskedValues(~isnan(maskedValues));

        if isempty(maskedValues)
            connectivityMatrix(i, j) = 0;
        else
            switch metric
                case 'mean'
                    connectivityMatrix(i, j) = mean(maskedValues, 'omitnan');
                case 'max'
                    connectivityMatrix(i, j) = max(maskedValues);
                case 'volume'
                    connectivityMatrix(i, j) = sum(maskedValues, 'omitnan') * voxelVolume;
            end
        end
    end
end

%% Apply normalization
if normalizeRows
    rowMax = max(connectivityMatrix, [], 2);
    rowMax(rowMax == 0) = 1; % Avoid division by zero
    connectivityMatrix = connectivityMatrix ./ rowMax;
end

if normalizeCols
    colMax = max(connectivityMatrix, [], 1);
    colMax(colMax == 0) = 1;
    connectivityMatrix = connectivityMatrix ./ colMax;
end

%% Create figure and plot heatmap
figHandle = figure('Color', 'w');

imagesc(connectivityMatrix);

% Set colormap
if ischar(cmapName)
    colormap(cmapName);
else
    colormap(cmapName);
end

% Add colorbar
cb = colorbar;
switch metric
    case 'mean'
        cb.Label.String = 'Mean projection density';
    case 'max'
        cb.Label.String = 'Max projection density';
    case 'volume'
        cb.Label.String = 'Total projection volume (mm³)';
end

% Set axis labels
set(gca, 'XTick', 1:nTargets, 'XTickLabel', targetRegions);
set(gca, 'YTick', 1:nSources, 'YTickLabel', sourceRegions);
xtickangle(45);

xlabel('Target Region');
ylabel('Source Region');

if ~isempty(plotTitle)
    title(plotTitle);
end

%% Annotate cells with values
if annotate
    % Normalize matrix for text color determination
    normMatrix = connectivityMatrix - min(connectivityMatrix(:));
    if max(connectivityMatrix(:)) > min(connectivityMatrix(:))
        normMatrix = normMatrix / (max(connectivityMatrix(:)) - min(connectivityMatrix(:)));
    end

    for i = 1:nSources
        for j = 1:nTargets
            value = connectivityMatrix(i, j);

            % Use white text on dark cells, black on light
            if normMatrix(i, j) < 0.5
                textColor = 'w';
            else
                textColor = 'k';
            end

            % Format based on magnitude
            if strcmp(metric, 'volume')
                if abs(value) < 0.01
                    txt = sprintf('%.2e', value);
                else
                    txt = sprintf('%.3f', value);
                end
            else
                if abs(value) < 0.001
                    txt = sprintf('%.2e', value);
                else
                    txt = sprintf('%.4f', value);
                end
            end

            text(j, i, txt, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', textColor, 'FontSize', 8);
        end
    end
end

axis image;

end
