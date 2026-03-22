function experimentIDs = findConnectivityExperiments(regions, mouseLine, primaryInjection)

%% sanitize inputs
if nargin < 2
    mouseLine = '';
end

if nargin < 3 || isempty(primaryInjection)
    primaryInjection = true;
end

experimentIDs = [];

for iRegion = 1:size(regions,2)

    %% build the query URL
    baseURL = 'http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure';

    % add region of interest
    baseURL = [baseURL, '[injection_structures$eq', regions{iRegion}, ']'];

    % mouse line
    if ~isempty(mouseLine)
        baseURL = [baseURL, '[transgenic_lines$eq', mouseLine, ']'];
    end

    % primary structure or not
    if primaryInjection
        primary = 'true';
    else
        primary = 'false';
    end
    fullURL = [baseURL, '[primary_structure_only$eq', primary, ']'];

    %Get data from Allen
    experimentPage = urlread(fullURL);
    try
        result = jsondecode(experimentPage);
    catch
        keyboard;
    end
    if ~result.success
        fprintf('Query failed!\n%s\nAt URL: %s\n\n', result.msg, fullURL);
    end

    %% get the experiment IDs
    %experimentIDs = [experimentIDs, ones(1, length(result.msg))];
    for iID = 1:length(result.msg)
        if isfield(result.msg, 'id')
            experimentIDs = [experimentIDs, result.msg(iID).id];
        else
            %fprintf('Query empty!\n');
        end
    end
    fprintf('Found %d experiments in %s \n', length(result.msg), regions{iRegion})
end
end