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
    disp(['Loading ' num2str(nExpIDs) ' experiments...']);
    progressBarHandle = waitbar(0, 'Loading experiments...');

    for iExpID = 1:nExpIDs
        waitbar(iExpID / totalExperiments, progressBarHandle, sprintf('Loading experiment %d of %d', iExpID, nExpIDs)); 
        
        % create dir if it doesn't exist
        saveDir = [saveLocation, filesep, num2str(experimentIDs(iExpID))];
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        % summary (structure.ionizes)
        summaryFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'injectionSummary_all.mat'];
        if ~exist(summaryFilePath)
            status = bsv.fetchConnectivitySummary(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
            if ~true
                continue;
            end
        end
         load(summaryFilePath,'injectionInfo'); % MC specified what to load
        
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
            if size(combinedInjectionInfo.(fieldNames{iFieldName}),1) == 1
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

        [sorted_values, sort_indices] = st.acronym(st.id==combinedInjectionInfo.structure_id);
        [groups, ~, groupID] = unique(sorted_values);
        groupID_original_order = zeros(size(groupID));
        groupID_original_order(sort_indices) = groupID;
    
    else
        warning('grouping method not recognized - skipping grouping. ')
        groups = ones(nExpIDs, 1);
        groupID_original_order = ones(nExpIDs, 1);
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
        waitbar(iExpID / nExpIDs, progressBarHandle, sprintf('Getting raw image %d of %d', iExpID, nExpIDs));

        thisGroup = groupID_original_order(iExpID);

        % raw data
        rawFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'density.raw'];
        if ~exist(rawFilePath)
            status = bsv.fetchConnectivityImages(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
            if ~status
                continue;
            end
        end

        fid = fopen(rawFilePath, 'r', 'l');
        experiment_projection = fread(fid, prod(projectionGridSize), 'float');
        fclose(fid);
        experiment_projection = reshape(experiment_projection, projectionGridSize);


        % Sum or average the projections
        % if strcmp(normalizationMethod, 'injectionIntensity')
        %     experiment_projection = experiment_projection / (combinedInjectionInfo.sum_projection_pixel_intensity(expID)...
        % /combinedInjectionInfo.sum_pixels(expID));
        % else
        %     experiment_projection = experiment_projection;
        % end

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
            individualProjections(:,:,:,iExpID) = experiment_projection;
        end
        combinedProjection(:,:,:,thisGroup) = combinedProjection(:,:,:,thisGroup) + experiment_projection;


    end
    close(progressBarHandle);
    
    % normalize images to number of projections - also need to add
    % normalization to injection intensity / volume
    for iGroup = 1:numberOfGroups
        theseGroups = numel(find(groupID_original_order==iGroup));
        combinedProjection(:,:,:,numberOfGroups) = combinedProjection(:,:,:,numberOfGroups) ./ theseGroups;
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