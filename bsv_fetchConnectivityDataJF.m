function [combinedProjection, combinedInjectionInfo] = bsv_fetchConnectivityDataJF(experimentIDs, saveLocation, fileName, ...
    normalizationMethod, subtractOtherHemisphere, groupingMethod)

if nargin < 6 || isempty(groupingMethod)% group experiments by brain region (groupingMethod = 'region') or not
    groupingMethod = 'NaN';
end

%% fetch data

projectionGridSize = [132, 80, 114];

% filepaths 
filePath_imgs = [saveLocation, filesep, fileName, '_', normalizationMethod, '_sub', num2str(subtractOtherHemisphere), '.mat'];
filePath_injectionSummary = [saveLocation, filesep, fileName, '_injectionSummary.mat'];

if ~exist(filePath_imgs, 'file') || isempty(fileName)
    combinedInjectionInfo = table;
    for iExpID = 1:size(experimentIDs, 2)
        
        % create dir if it doesn't exist
        saveDir = [saveLocation, filesep, num2str(experimentIDs(iExpID))];
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        % summary (structure.ionizes)
        summaryFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'injectionSummary.mat'];
        if ~exist(summaryFilePath)
            status = bsv_fetchConnectivitySummary(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
            if ~true
                continue;
            end
        end
        load(summaryFilePath)
        
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
            combinedInjectionInfo.(fieldNames{iFieldName}) = [combinedInjectionInfo.(fieldNames{iFieldName}); injectionInfo.(fieldNames{iFieldName})];
        end

    end
    
    % grouping method 
    if strcmp(groupingMethod, 'AP')
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
        groups = ones(size(experimentIDs,2), 1);
        groupID_original_order = ones(size(experimentIDs,2), 1);
    else
        warning('grouping method not recognized - skipping grouping. ')
        groups = ones(size(experimentIDs,2), 1);
        groupID_original_order = ones(size(experimentIDs,2), 1);
    end
    

    numberOfGroups = numel(unique(groups));
    combinedProjection = zeros([projectionGridSize, numberOfGroups]);
    
    % get raw images 
    for iExpID = 1:size(experimentIDs, 2)

        thisGroup = groupID_original_order(iExpID);

        % raw data
        rawFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'density.raw'];
        if ~exist(rawFilePath)
            status = bsv_fetchConnectivityImages(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
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
            % debug
            %figure(); imagesc(squeeze(experiment_projection(88,:,:)))
            %figure(); imagesc(squeeze(experiment_projection_tmp(88,:,:)))
            %figure(); imagesc(squeeze(experiment_projection(88,:,:))-squeeze(experiment_projection_tmp(88,:,:)))

            experiment_projection = experiment_projection_tmp;
        end

        combinedProjection(:,:,:,thisGroup) = combinedProjection(:,:,:,thisGroup) + experiment_projection;


    end
    
    % normalize images to number of projections 
    for iGroup = 1:numberOfGroups
        theseGroups = numel(find(groupID_original_order==iGroup));
        combinedProjection(:,:,:,numberOfGroups) = combinedProjection(:,:,:,numberOfGroups) ./ theseGroups;
    end

    % save files 
    if ~isempty(fileName)
        save(filePath_imgs, 'combinedProjection')
        save(filePath_injectionSummary, 'combinedInjectionInfo')
    end
else
    load(filePath_imgs)
    readtable(filePath_injectionSummary)
end

end