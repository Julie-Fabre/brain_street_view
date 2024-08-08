function [combinedProjection, combinedInjectionInfo, individualProjections] = fetchConnectivityData(experimentIDs, saveLocation, fileName, ...
    normalizationMethod, subtractOtherHemisphere, groupingMethod, allenAtlasPath, loadAll)

if nargin < 6 || isempty(groupingMethod) || nargin < 7 || isempty(allenAtlasPath) % group experiments by brain region (groupingMethod = 'region') or not
    groupingMethod = 'NaN';
else
    st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
end

%% fetch data

projectionGridSize = [132, 80, 114];

% filepaths
filePath_imgs = [saveLocation, filesep, fileName, '_', normalizationMethod, '_sub', num2str(subtractOtherHemisphere), '.mat'];
filePath_injectionSummary = [saveLocation, filesep, fileName, '_injectionSummary.mat'];

nExpIDs = size(experimentIDs, 2); % added by Matteo

if ~exist(filePath_imgs, 'file') || isempty(fileName)
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
            mkdir(saveDir);
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

    end
    close(progressBarHandle);

    % grouping method
    if strcmp(groupingMethod, 'brainRegion')
        % outputAcronyms = arrayfun(@(id) st.acronym{st.id == id}, combinedInjectionInfo.structure_id, 'UniformOutput', false); % Commented out by Matteo because it's unused
        [sorted_values, sort_indices] = st.acronym(st.id == combinedInjectionInfo.structure_id);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'AP')
        [sorted_values, sort_indices] = sort(combinedInjectionInfo.max_voxel_x);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'ML')
        [sorted_values, sort_indices] = sort(combinedInjectionInfo.max_voxel_z);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'DV')
        [sorted_values, sort_indices] = sort(combinedInjectionInfo.max_voxel_y);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    elseif strcmp(groupingMethod, 'NaN') || isempty(groupingMethod)
        groups = ones(size(experimentIDs, 2), 1);
        groupID_original_order = ones(size(experimentIDs, 2), 1);
    else
        warning('grouping method not recognized - skipping grouping. ')
        groups = ones(size(experimentIDs, 2), 1);
        groupID_original_order = ones(size(experimentIDs, 2), 1);
    end


    numberOfGroups = numel(unique(groups));
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
        rawFilePath = fullfile(saveLocation, filesep, num2str(currExpID), filesep, 'density.raw');
        if ~exist(rawFilePath, 'file')
            % fetch and save raw data if it isn't already on disk
            status = bsv.fetchConnectivityImages(currExpID, [saveLocation, filesep, num2str(currExpID)]);
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
        if strcmp(normalizationMethod, 'injectionVolume')

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
            
            % Extract projection volumes for hemisphere 3 and hemispheres 1/2
            vol3 = combinedInjectionInfo.projection_volume(hem3Indices);
            vol12 = combinedInjectionInfo.projection_volume(hem12Indices);
            
            % Check conditions for valid calculation
            if (numel(vol3) == 1 && numel(vol12) == 1) || ...
               (numel(unique(vol3)) == 1 && numel(unique(vol12)) == 1)
                % there should be one unique entry for each. sometimes the values are duplicated. no diea why but this is fine
                
                % If values are duplicated, use the first one
                vol3 = vol3(1);
                vol12 = vol12(1);
                
                % Calculate volume difference and adjusted volume
                volumeDifference = vol12 - vol3;
                adjustedVolume = vol12 + volumeDifference;
            else
                warning('Unexpected number of entries for structure 997 in experiment %d. Expected one entry for hemisphere 3 and one for hemisphere 1/2.', currExpID);
                keyboard;
                adjustedVolume = NaN;
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
        combinedProjection(:, :, :, numberOfGroups) = combinedProjection(:, :, :, numberOfGroups) ./ theseGroups;
    end

    % save files
    if ~isempty(fileName)
        disp('Saving combined results...');
        save(filePath_imgs, 'combinedProjection')
        save(filePath_injectionSummary, 'combinedInjectionInfo')
    end
else
    disp('Loading existing data...');
    load(filePath_imgs, 'combinedProjection')
    readtable(filePath_injectionSummary)
end
disp('Data processing complete.');
end