function [combinedProjection, combinedInjectionInfo, individualProjections, experimentRegionInfo] = fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    normalizationMethod, subtractOtherHemisphere, groupingMethod, allenAtlasPath, loadAll, inputRegions, regionGroups, exportMetadata, reload)

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
    exportMetadata = false;
end
if nargin < 12 || isempty(reload)
    reload = false;
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

% filepaths
filePath_imgs = [saveLocation, filesep, fileName, '_', normalizationMethod, '_sub', num2str(subtractOtherHemisphere), '.mat'];
filePath_injectionSummary = [saveLocation, filesep, fileName, '_injectionSummary.mat'];

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


if ~exist(filePath_imgs, 'file') || isempty(fileName) || reload
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
        for iExpID = 1:nExpIDs
            regionAbbrev = primaryStructure_abbreviation{iExpID};
            
            % Find which input region this experiment corresponds to
            regionIdx = 0;
            for iRegion = 1:length(inputRegions)
                if strcmp(regionAbbrev, inputRegions{iRegion}) || contains(regionAbbrev, inputRegions{iRegion})
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

            % Find indices for the current experiment
            expIndices = combinedInjectionInfo.experimentID == currExpID;

            % Find indices for structure 997 (root-> so it will contain all pther structures) within this experiment
            % another option is to just use volume for injection structure of
            % interest. The best solution will depend on the application. 
            structureIndices = expIndices & combinedInjectionInfo.structure_id == 997;
            
            % Calculate the difference between hemisphere 3 (both hemipsheres: contains summed 
            % value for structure in hemisphere 1 + 2 and hemispheres 1/2
            hem3Indices = structureIndices & [combinedInjectionInfo.hemisphere_id] == 3;
            hem12Indices = structureIndices & ([combinedInjectionInfo.hemisphere_id] == 1 | [combinedInjectionInfo.hemisphere_id] == 2);
            
          if strcmp(normalizationMethod, 'injectionVolume')
            % Extract projection volumes for hemisphere 3 and hemispheres 1/2
            vol3 = combinedInjectionInfo.projection_volume(hem3Indices);
            vol12 = combinedInjectionInfo.projection_volume(hem12Indices);
          else
            vol3 = combinedInjectionInfo.sum_pixel_intensity(hem3Indices);
            vol12 = combinedInjectionInfo.sum_pixel_intensity(hem12Indices);
          end
            
            % Handle multiple entries by taking appropriate values
            % vol3 should be a single value (hemisphere 3 = both hemispheres)
            % vol12 can have multiple entries (separate left/right hemisphere values)
            
            vol3_single = max(vol3);  % Use max in case of duplicates
            vol12_sum = sum(vol12);   % Sum the hemisphere values
            
            % Check if data is valid (non-zero and reasonable)
            if vol3_single > 0 && vol12_sum > 0
                % Calculate volume difference and adjusted volume
                volumeDifference = vol12_sum - vol3_single;
                adjustedVolume = vol12_sum + volumeDifference;
            else
                fprintf('Warning: Invalid volume data for experiment %d (vol3=%.6f, vol12_sum=%.6f), using fallback\n', currExpID, vol3_single, vol12_sum);
                adjustedVolume = max(vol3_single, vol12_sum);
                if adjustedVolume == 0
                    adjustedVolume = 1; % Prevent division by zero
                end
            end
                        
            % normalize to adjusted volume
            experiment_projection = experiment_projection / adjustedVolume;
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
    
    % save files
    if ~isempty(fileName)
        disp('Saving combined results...');
        save(filePath_imgs, 'combinedProjection')
        save(filePath_injectionSummary, 'combinedInjectionInfo')
        save([saveLocation, filesep, fileName, '_experimentRegionInfo.mat'], 'experimentRegionInfo')
    end
else
    disp('Loading existing data...');
    load(filePath_imgs, 'combinedProjection')
    combinedInjectionInfo = load(filePath_injectionSummary);
    combinedInjectionInfo = combinedInjectionInfo.combinedInjectionInfo;
    
    % Initialize individualProjections when loading from cache
    if loadAll
        % When loading from cache, we don't have individual projections stored
        % so we'll initialize with empty array
        individualProjections = [];
        warning('Individual projections not available when loading from cache. Re-run without cache to get individual projections.');
    else
        individualProjections = '';
    end
    
    % Load experiment region info if it exists
    experimentRegionInfoPath = [saveLocation, filesep, fileName, '_experimentRegionInfo.mat'];
    if exist(experimentRegionInfoPath, 'file')
        load(experimentRegionInfoPath, 'experimentRegionInfo')
    else
        % Create from Allen atlas info if cached data doesn't have it
        experimentRegionInfo.experimentIDs = experimentIDs;
        nExpIDs = length(experimentIDs);
        experimentRegionInfo.primaryStructure_ID = zeros(nExpIDs, 1);
        experimentRegionInfo.primaryStructure_abbreviation = cell(nExpIDs, 1);
        
        % Load allen atlas info
        currentFileLocation =  which(mfilename('fullpath'));
        currentFileLocation_Pa = fileparts(fileparts(currentFileLocation));
        allenAtlasProjection_info = readtable(fullfile(currentFileLocation_Pa, 'docs', 'allenAtlasProjection_info.csv'), 'VariableNamingRule','modify');
        
        for iExpID = 1:nExpIDs
            currExpID = experimentIDs(iExpID);
            experimentRegionInfo.primaryStructure_ID(iExpID) = allenAtlasProjection_info.structure_id(allenAtlasProjection_info.id == currExpID);
            experimentRegionInfo.primaryStructure_abbreviation{iExpID} = allenAtlasProjection_info.structure_abbrev{allenAtlasProjection_info.id == currExpID};
        end
    end
    
    % Export metadata to CSV if requested (for loaded data)
    if exportMetadata
        % Use default filename if none provided
        csvFileName_base = fileName;
        if isempty(csvFileName_base) || strcmp(csvFileName_base, '')
            csvFileName_base = 'connectivity_data';
            fprintf('No filename provided for metadata export, using default: %s\n', csvFileName_base);
        end
        
        csvFileName = [saveLocation, filesep, csvFileName_base, '_metadata.csv'];
        
        % Check if CSV already exists
        if exist(csvFileName, 'file')
            fprintf('Metadata CSV already exists: %s\n', csvFileName);
        else
            fprintf('Exporting experiment metadata to CSV (from loaded data)...\n');
            
            % Create metadata table from loaded data
            metadataTable = table();
            metadataTable.ExperimentID = experimentIDs(:);
            
            % Get experiment info from Allen Atlas data
            expInfo = allenAtlasProjection_info(ismember(allenAtlasProjection_info.id, experimentIDs), :);
            [~, idx] = ismember(experimentIDs, expInfo.id);
            expInfo = expInfo(idx, :);
            
            % Mouse line/genotype
            metadataTable.MouseLine = expInfo.transgenic_line;
            emptyIdx = cellfun(@isempty, metadataTable.MouseLine) | strcmp(metadataTable.MouseLine, '""') | strcmp(metadataTable.MouseLine, '');
            metadataTable.MouseLine(emptyIdx) = {'Wild-type'};
            
            % Injection region info
            if isfield(experimentRegionInfo, 'primaryStructure_abbreviation')
                metadataTable.InjectionRegion = experimentRegionInfo.primaryStructure_abbreviation(:);
            else
                metadataTable.InjectionRegion = repmat({'Unknown'}, length(experimentIDs), 1);
            end
            
            if isfield(experimentRegionInfo, 'primaryStructure_ID')
                metadataTable.InjectionRegionID = experimentRegionInfo.primaryStructure_ID(:);
            else
                metadataTable.InjectionRegionID = zeros(length(experimentIDs), 1);
            end
            
            % Injection coordinates (if available)
            if isfield(experimentRegionInfo, 'primaryStructure_AP')
                metadataTable.InjectionAP = experimentRegionInfo.primaryStructure_AP(:);
                metadataTable.InjectionDV = experimentRegionInfo.primaryStructure_DV(:);
                metadataTable.InjectionML = experimentRegionInfo.primaryStructure_ML(:);
            else
                metadataTable.InjectionAP = zeros(length(experimentIDs), 1);
                metadataTable.InjectionDV = zeros(length(experimentIDs), 1);
                metadataTable.InjectionML = zeros(length(experimentIDs), 1);
            end
            
            % Get injection volume and intensity from loaded combinedInjectionInfo
            injectionVolumes = zeros(length(experimentIDs), 1);
            injectionIntensities = zeros(length(experimentIDs), 1);
            
            for iExp = 1:length(experimentIDs)
                expID = experimentIDs(iExp);
                
                % Find injection data for this experiment
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
            
            % Sex and age if available
            if ismember('sex', expInfo.Properties.VariableNames)
                metadataTable.Sex = expInfo.sex;
            end
            if ismember('age', expInfo.Properties.VariableNames)
                metadataTable.Age = expInfo.age;
            end
            
            % Save CSV file
            writetable(metadataTable, csvFileName);
            fprintf('Metadata exported to: %s\n', csvFileName);
        end
    end
end
disp('Data processing complete.');
end