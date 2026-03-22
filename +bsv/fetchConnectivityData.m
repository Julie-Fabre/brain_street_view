function [combinedProjection, combinedInjectionInfo, individualProjections, experimentRegionInfo] = fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    normalizationMethod, subtractOtherHemisphere, groupingMethod, allenAtlasPath, loadAll, inputRegions, regionGroups, exportMetadata, reload, atlasType, atlasResolution)

if nargin < 6 || isempty(groupingMethod) || nargin < 7 || isempty(allenAtlasPath) % group experiments by brain region (groupingMethod = 'region') or not
    groupingMethod = 'NaN';
%else
    % st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
    % I will need this to implement some selection of regions by hierarchy.
    % 
end

% Handle optional input regions and region groups
if nargin < 9 || isempty(inputRegions)
    inputRegions = [];
end
if nargin < 10 || isempty(regionGroups)
    regionGroups = [];
end
if nargin < 11 || isempty(exportMetadata)
    exportMetadata = true;
end
if nargin < 12 || isempty(reload)
    reload = false;
end
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

%% fetch data

projectionGridSize = [132, 80, 114];

% Create main save directory if it doesn't exist
if ~exist(saveLocation, 'dir')
    try
        mkdir(saveLocation);
    catch ME
        warning('Could not create directory %s: %s\nUsing current directory instead.', saveLocation, ME.message);
        saveLocation = pwd;
    end
end

% filepaths - include atlas info in filename
atlasStr = sprintf('_%s_%dum', atlasType, atlasResolution);
filePath_imgs = [saveLocation, filesep, fileName, '_', normalizationMethod, '_sub', num2str(subtractOtherHemisphere), atlasStr, '.mat'];
filePath_injectionSummary = [saveLocation, filesep, fileName, '_injectionSummary', atlasStr, '.mat'];

% load summary from the Allen
currentFileLocation =  which(mfilename('fullpath'));
currentFileLocation_Pa = fileparts(fileparts(currentFileLocation)); % twp directories up from current file location
allenAtlasProjection_info = readtable(fullfile(currentFileLocation_Pa, 'docs', 'allenAtlasProjection_info.csv'), 'VariableNamingRule','modify');

% intialize arrays
nExpIDs = size(experimentIDs, 2); % added by Matteo
primaryStructure_AP = zeros(nExpIDs, 1);
primaryStructure_DV = zeros(nExpIDs, 1);
primaryStructure_ML = zeros(nExpIDs, 1);
primaryStructure_ID = zeros(nExpIDs, 1);
primaryStructure_abbreviation = cell(nExpIDs, 1);



    combinedInjectionInfo = table;
    
    % display progress
    disp(['Loading ', num2str(nExpIDs), ' experiments...']);
    progressBarHandle = waitbar(0, 'Loading experiments...');

    for iExpID = 1:nExpIDs
        waitbar(iExpID/nExpIDs, progressBarHandle, sprintf('Loading experiment %d of %d', iExpID, nExpIDs));

        currExpID = experimentIDs(iExpID);
        % create dir if it doesn't exist
        saveDir = [saveLocation, filesep, num2str(currExpID)];
        if ~exist(saveDir, 'dir')
            try
                [status, msg] = mkdir(saveDir);
                if ~status
                    error('Failed to create directory %s: %s', saveDir, msg);
                end
            catch ME
                warning('Could not create experiment directory %s: %s\nSkipping download for experiment %d.', saveDir, ME.message, currExpID);
                continue;
            end
        end

        % summary (structure.ionizes)
        summaryFilePath = fullfile(saveLocation, num2str(currExpID), 'injectionSummary_all.mat');
        if ~exist(summaryFilePath, 'file')
            status = bsv.fetchConnectivitySummary(currExpID, [saveLocation, filesep, num2str(currExpID)]);
            if ~status
                continue;
            end
        end
        load(summaryFilePath, 'injectionInfo'); % MC specified what to load

        % put all strcuture.ionizes summary into one table
        % Initialize combinedInjectionInfo if it doesn't exist
        if iExpID == 1
            combinedInjectionInfo = struct();
            fieldNames = fieldnames(injectionInfo);
        end

        for iFieldName = 1:numel(fieldNames)
            % Ensure the field exists in combinedInjectionInfo
            if ~isfield(combinedInjectionInfo, fieldNames{iFieldName})
                combinedInjectionInfo.(fieldNames{iFieldName}) = [];
            end

            % Concatenate the data
            if size(combinedInjectionInfo.(fieldNames{iFieldName}), 1) == 1
                combinedInjectionInfo.(fieldNames{iFieldName}) = [combinedInjectionInfo.(fieldNames{iFieldName}), [injectionInfo.(fieldNames{iFieldName})]];
            else
                combinedInjectionInfo.(fieldNames{iFieldName}) = [combinedInjectionInfo.(fieldNames{iFieldName}); [injectionInfo.(fieldNames{iFieldName})]];
            end
        end
        currInjectionInfo = struct2table(injectionInfo);
        primaryStructure_ID(iExpID) = allenAtlasProjection_info.structure_id(allenAtlasProjection_info.id == currExpID); % this should hold true 
        % in most cases. the best way would be to take into account group hierarchies but they
        % make no sense to me for instance VISp is at level 7, CP and GPe at level 6, SNr at level
        % 5 but conceptually they are at the same level to me...
        primaryStructure_abbreviation{iExpID} = allenAtlasProjection_info.structure_abbrev{allenAtlasProjection_info.id == currExpID};
        currInjectionCoordinates = str2num(allenAtlasProjection_info.injection_coordinates{allenAtlasProjection_info.id == currExpID});
        primaryStructure_AP(iExpID) = currInjectionCoordinates(1);
        primaryStructure_DV(iExpID) = currInjectionCoordinates(2);
        primaryStructure_ML(iExpID) = currInjectionCoordinates(3);
    end
    close(progressBarHandle);

    %% Display experiment summary
    % Filter for current experiments
    expInfo = allenAtlasProjection_info(ismember(allenAtlasProjection_info.id, experimentIDs), :);
    
    fprintf('\n📊 EXPERIMENT SUMMARY\n');
    fprintf('Total experiments loaded: %d\n', length(experimentIDs));
    
    % Mouse genotype distribution
    fprintf('\n🧬 Mouse genotype distribution:\n');
    genotypes = expInfo.transgenic_line;
    
    % Handle empty genotypes
    genotypes(cellfun(@isempty, genotypes)) = {'Wild-type'};
    genotypes(strcmp(genotypes, '""') | strcmp(genotypes, '')) = {'Wild-type'};
    
    [uniqueGenotypes, ~, genotypeIdx] = unique(genotypes);
    for i = 1:length(uniqueGenotypes)
        count = sum(genotypeIdx == i);
        fprintf('  • %s: %d experiments\n', uniqueGenotypes{i}, count);
    end
    
    % Brain region distribution
    fprintf('\n🧠 Brain region distribution:\n');
    regions = expInfo.structure_abbrev;
    [uniqueRegions, ~, regionIdx] = unique(regions);
    for i = 1:length(uniqueRegions)
        count = sum(regionIdx == i);
        fprintf('  • %s: %d experiments\n', uniqueRegions{i}, count);
    end
    fprintf('─────────────────────────────────\n\n');

    % % grouping method 
    if strcmp(groupingMethod, 'brainRegion')
        % because of hierarchy differences, this doesn't make much sense
        % currently, not implemented 
        warning('grouping method not implemented yet - skipping grouping. ')
        % [sorted_values, sort_indices] = st.acronym(st.id == combinedInjectionInfo.structure_id);
        % [groups, ~, groupID] = unique(primaryStructure_ID);
        % groupID_original_order = zeros(size(groupID));
        % groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'AP') % group by AP value
        [sorted_values, sort_indices] = sort(primaryStructure_AP);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'ML') % group by ML value
        [sorted_values, sort_indices] = sort(primaryStructure_ML);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'DV') % group by DV value
        [sorted_values, sort_indices] = sort(primaryStructure_DV);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'NaN') || isempty(groupingMethod) % no grouping
        groups = ones(size(experimentIDs, 2), 1);
        groupID_original_order = ones(size(experimentIDs, 2), 1);
    else % error
        warning('grouping method not recognized - skipping grouping. ')
        groups = ones(size(experimentIDs, 2), 1);
        groupID_original_order = ones(size(experimentIDs, 2), 1);
    end


    numberOfGroups = numel(unique(groups));
    
    % Override grouping if region-based grouping is requested
    if ~isempty(inputRegions) && ~isempty(regionGroups)
        fprintf('Setting up region-based grouping for %d input regions into %d groups\n', length(inputRegions), length(unique(regionGroups)));
        
        % Create region-to-group mapping
        uniqueRegionGroups = unique(regionGroups);
        nRegionGroups = length(uniqueRegionGroups);
        regionGroupsCellArray = cell(nRegionGroups, 1);
        
        for i = 1:nRegionGroups
            regionGroupsCellArray{i} = find(regionGroups == uniqueRegionGroups(i));
        end
        
        % Map each experiment to its region group
        regionBasedGroupID = zeros(nExpIDs, 1);
        unmatchedExperiments = {};
        for iExpID = 1:nExpIDs
            regionAbbrev = primaryStructure_abbreviation{iExpID};
            
            % Find which input region this experiment corresponds to
            regionIdx = 0;
            for iRegion = 1:length(inputRegions)
                if strcmp(regionAbbrev, inputRegions{iRegion}) || startsWith(regionAbbrev, inputRegions{iRegion})
                    regionIdx = iRegion;
                    break;
                end
            end
            
            % Find which region group this region belongs to
            if regionIdx > 0
                for iRegionGroup = 1:nRegionGroups
                    if any(regionGroupsCellArray{iRegionGroup} == regionIdx)
                        regionBasedGroupID(iExpID) = iRegionGroup;
                        break;
                    end
                end
            else
                % Log unmatched experiments
                unmatchedExperiments{end+1} = sprintf('Exp %d: %s', experimentIDs(iExpID), regionAbbrev);
            end
        end
        
        % Report unmatched experiments
        if ~isempty(unmatchedExperiments)
            fprintf('\nWarning: %d experiments did not match any input region:\n', length(unmatchedExperiments));
            for i = 1:min(10, length(unmatchedExperiments))
                fprintf('  %s\n', unmatchedExperiments{i});
            end
            if length(unmatchedExperiments) > 10
                fprintf('  ... and %d more\n', length(unmatchedExperiments) - 10);
            end
        end
        
        % Override the grouping variables
        numberOfGroups = nRegionGroups;
        groupID_original_order = regionBasedGroupID;
        
        fprintf('Region-based grouping: %d experiments mapped to %d groups\n', sum(regionBasedGroupID > 0), numberOfGroups);
    end
    
    combinedProjection = zeros([projectionGridSize, numberOfGroups]);

    if loadAll
        individualProjections = zeros([projectionGridSize, size(experimentIDs, 2)]);
    else
        individualProjections = '';
    end

    % get raw images
    disp('Getting raw images...');
    progressBarHandle = waitbar(0, 'Getting raw images...');
    for iExpID = 1:nExpIDs
        waitbar(iExpID/nExpIDs, progressBarHandle, sprintf('Getting raw image %d of %d', iExpID, nExpIDs));

        currGroup = groupID_original_order(iExpID);
        currExpID = experimentIDs(iExpID);
        
        % Skip experiments that don't belong to any group
        if currGroup == 0
            continue;
        end

        % raw data
        rawFilePath = fullfile(saveLocation, num2str(currExpID), 'density.raw');
        if ~exist(rawFilePath, 'file')
            % fetch and save raw data if it isn't already on disk
            status = bsv.fetchConnectivityImages(currExpID, fullfile(saveLocation, num2str(currExpID)));
            if ~status
                continue;
            end
        end
        
        % load raw data file
        fid = fopen(rawFilePath, 'r', 'l');
        experiment_projection = fread(fid, prod(projectionGridSize), 'float');
        fclose(fid);
        experiment_projection = reshape(experiment_projection, projectionGridSize);
        
        % normalize projection values
        if strcmp(normalizationMethod, 'injectionVolume') || strcmp(normalizationMethod, 'injectionIntensity')

            % Find indices for the current experiment and structure 997 (root)
            expIndices = combinedInjectionInfo.experimentID == currExpID;
            structureIndices = expIndices & combinedInjectionInfo.structure_id == 997;

            % Use hemisphere 3 (both hemispheres combined)
            hem3Indices = structureIndices & [combinedInjectionInfo.hemisphere_id] == 3;

            % Both methods normalize by projection_volume (mm^3).
            % The density.raw values are already fractional (0-1), so dividing
            % by injection volume makes experiments with different injection
            % sizes comparable.
            injVol = max(combinedInjectionInfo.projection_volume(hem3Indices));

            if isempty(injVol) || injVol == 0
                fprintf('Warning: Invalid injection volume for experiment %d, skipping normalization\n', currExpID);
                injVol = 1;
            end

            experiment_projection = experiment_projection / injVol;
        end

        if subtractOtherHemisphere
            experiment_projection_tmp = zeros(132, 80, 114);
            if injectionInfo.max_voxel_z <= 114 / 2 %left
                for iML = 1:57
                    experiment_projection_tmp(:, :, iML) = [experiment_projection(:, :, iML) - experiment_projection(:, :, 114-iML+1)];
                end
            elseif injectionInfo.max_voxel_z >= 114 / 2 %right
                for iML = 1:57
                    experiment_projection_tmp(:, :, 114-iML+1) = [experiment_projection(:, :, 114-iML+1) - experiment_projection(:, :, iML)];
                end
            end
            experiment_projection = experiment_projection_tmp;
        end
        if loadAll
            individualProjections(:, :, :, iExpID) = experiment_projection;
        end
        combinedProjection(:, :, :, currGroup) = combinedProjection(:, :, :, currGroup) + experiment_projection;


    end
    close(progressBarHandle);

    % normalize images to number of projections - also need to add
    % normalization to injection intensity / volume
    for iGroup = 1:numberOfGroups
        theseGroups = numel(find(groupID_original_order == iGroup));
        combinedProjection(:, :, :, iGroup) = combinedProjection(:, :, :, iGroup) ./ theseGroups;
    end

    % Create experiment region info structure
    experimentRegionInfo.experimentIDs = experimentIDs;
    experimentRegionInfo.primaryStructure_ID = primaryStructure_ID;
    experimentRegionInfo.primaryStructure_abbreviation = primaryStructure_abbreviation;
    experimentRegionInfo.primaryStructure_AP = primaryStructure_AP;
    experimentRegionInfo.primaryStructure_DV = primaryStructure_DV;
    experimentRegionInfo.primaryStructure_ML = primaryStructure_ML;
    
    % Export metadata to CSV if requested
    if exportMetadata
        % Use default filename if none provided
        csvFileName_base = fileName;
        if isempty(csvFileName_base) || strcmp(csvFileName_base, '')
            csvFileName_base = 'connectivity_data';
            fprintf('No filename provided for metadata export, using default: %s\n', csvFileName_base);
        end
        fprintf('Exporting experiment metadata to CSV...\n');
        
        % Create metadata table
        metadataTable = table();
        
        % Basic experiment info
        metadataTable.ExperimentID = experimentIDs(:);
        
        % Get experiment info from Allen Atlas data
        expInfo = allenAtlasProjection_info(ismember(allenAtlasProjection_info.id, experimentIDs), :);
        
        % Ensure the order matches experimentIDs
        [~, idx] = ismember(experimentIDs, expInfo.id);
        expInfo = expInfo(idx, :);
        
        % Mouse line/genotype
        metadataTable.MouseLine = expInfo.transgenic_line;
        % Clean up empty genotypes
        emptyIdx = cellfun(@isempty, metadataTable.MouseLine) | strcmp(metadataTable.MouseLine, '""') | strcmp(metadataTable.MouseLine, '');
        metadataTable.MouseLine(emptyIdx) = {'Wild-type'};
        
        % Injection region
        metadataTable.InjectionRegion = primaryStructure_abbreviation(:);
        metadataTable.InjectionRegionID = primaryStructure_ID(:);
        
        % Injection coordinates
        metadataTable.InjectionAP = primaryStructure_AP(:);
        metadataTable.InjectionDV = primaryStructure_DV(:);
        metadataTable.InjectionML = primaryStructure_ML(:);
        
        % Get injection volume and intensity from combinedInjectionInfo
        % We need to extract the primary injection site data
        injectionVolumes = zeros(length(experimentIDs), 1);
        injectionIntensities = zeros(length(experimentIDs), 1);
        
        for iExp = 1:length(experimentIDs)
            expID = experimentIDs(iExp);
            
            % Find injection data for this experiment (structure_id = 997, hemisphere_id = 3)
            expIdx = combinedInjectionInfo.experimentID == expID;
            structIdx = combinedInjectionInfo.structure_id == 997;
            hemIdx = combinedInjectionInfo.hemisphere_id == 3;
            injIdx = expIdx & structIdx & hemIdx;
            
            if any(injIdx)
                % Handle multiple entries by taking the first/max value
                volumes = combinedInjectionInfo.projection_volume(injIdx);
                intensities = combinedInjectionInfo.sum_pixel_intensity(injIdx);
                
                if length(volumes) > 1
                    fprintf('Warning: Multiple injection entries found for experiment %d, using max values\n', expID);
                    injectionVolumes(iExp) = max(volumes);
                    injectionIntensities(iExp) = max(intensities);
                else
                    injectionVolumes(iExp) = volumes;
                    injectionIntensities(iExp) = intensities;
                end
            end
        end
        
        metadataTable.InjectionVolume = injectionVolumes;
        metadataTable.InjectionIntensity = injectionIntensities;
        
        % Sex information if available
        if ismember('sex', expInfo.Properties.VariableNames)
            metadataTable.Sex = expInfo.sex;
        end
        
        % Age information if available
        if ismember('age', expInfo.Properties.VariableNames)
            metadataTable.Age = expInfo.age;
        end
        
        % Save CSV file
        csvFileName = [saveLocation, filesep, csvFileName_base, '_metadata.csv'];
        writetable(metadataTable, csvFileName);
        fprintf('Metadata exported to: %s\n', csvFileName);
    end
disp('Data processing complete.');
end