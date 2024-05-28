function [combinedProjection, combinedInjectionInfo] = nsv_fetchConnectivityData(experimentIDs, saveLocation, fileName, normalizationMethod, subtractOtherHemisphere)

%% fetch data

projectionGridSize = [132, 80, 114];

combinedProjection = zeros(projectionGridSize);
combinedInjectionInfo = table;

filePath_imgs = [saveLocation, filesep, fileName, '_', normalizationMethod, '_sub', num2str(subtractOtherHemisphere), '.mat'];
filePath_injectionSummary = [saveLocation, filesep, fileName, '_injectionSummary.mat'];

if ~exist(filePath_imgs, 'file') || isempty(fileName)
    for iExpID = 1:size(experimentIDs, 2)

        mkdir([saveLocation, filesep, num2str(experimentIDs(iExpID))]);

        % raw data
        rawFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'density.raw'];
        if ~exist(rawFilePath)
            status = nsv_fetchConnectivityImages(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
            if ~status
                continue;
            end
        end

        fid = fopen(rawFilePath, 'r', 'l');
        experiment_projection = fread(fid, prod(projectionGridSize), 'float');
        fclose(fid);
        experiment_projection = reshape(experiment_projection, projectionGridSize);

        % summary (structure.ionizes
        summaryFilePath = [saveLocation, filesep, num2str(experimentIDs(iExpID)), filesep, 'injectionSummary.mat'];
        if ~exist(summaryFilePath)
            status = nsv_fetchConnectivitySummary(experimentIDs(iExpID), [saveLocation, filesep, num2str(experimentIDs(iExpID))]);
            if ~true
                continue;
            end
        end
        load(summaryFilePath)


        % Sum or average the projections
        % if strcmp(normalizationMethod, 'injectionIntensity')
        %     experiment_projection = experiment_projection / (injectionInfo.sum_projection_pixel_intensity/injectionInfo.sum_pixels);
        % else
        %     experiment_projection = experiment_projection;
        % end
        % 
        % if subtractOtherHemisphere
        %     if injectionInfo.hemisphere_id == 1 %left
        %         experiment_projection_tmp = zeros(132, 80, 114);
        %         for iML = 1:57
        %             experiment_projection_tmp(:, :, iML) = [experiment_projection(:, :, iML) - experiment_projection(:, :, 114-iML+1)];
        %         end
        %     elseif injectionInfo.hemisphere_id == 2 %right
        %         experiment_projection_tmp = zeros(132, 80, 114);
        %         for iML = 1:57
        %             experiment_projection_tmp(:, :, 114-iML+1) = [experiment_projection(:, :, 114-iML+1) - experiment_projection(:, :, iML)];
        %         end
        %     elseif injectionInfo.hemisphere_id == 3 %both
        %         % do nothing ?
        %     end
        % 
        % end

        combinedProjection = combinedProjection + experiment_projection;
        
        fieldNames = fieldnames(injectionInfo);
        for iFieldName = 1:size(fieldNames,1)
            combinedInjectionInfo.(fieldNames{iFieldName})(iExpID) = injectionInfo.(fieldNames{iFieldName});
        end

    end
    combinedProjection = combinedProjection ./ size(experimentIDs, 2);
    if ~isempty(fileName)
        save(filePath_imgs, 'combinedProjection')
        save(filePath_injectionSummary, 'combinedInjectionInfo')
    end
else
    load(filePath_imgs)
    readtable(filePath_injectionSummary)
end

end