function status = fetchConnectivitySummary(experimentID, saveFilePath)
% returns a MATLAB structure with one entry per injection structure and the
% following fields: 
% - injectionInfo.experimentID : the experiment ID.
% - structure_id : ID of the structure (e.g. 315 for Isocortex, 672 for caudoputamen).
% - hemisphere_id : ID of the hemisphere (1 = left, 2 = right, 3 = both)
%       if the entry is an injection & the hemisphere == 3, then the
%       injection itself spans both hemispheres:
%       https://community.brain-map.org/t/selecting-connectivity-data-with-injection-into-particular-hemisphere/2095 
% - is_injection : If true, numbers only include voxels from injection site.
%       If false, numbers only include voxels outside of the injection site.
%       When there are several entries for one experiment (=injection) with is_injection = true 
%       it means that the injection site spans multiple regions.
% - sum_pixels : Number of valid pixels in the structure in this experiment.
%       Valid pixels are those not manually annotated as invalid data.
% - sum_projection_pixels : Number of pixels identified as projecting in this
%       structure.
% - sum_pixel_intensity : Sum of intensity values in this structure.
% - sum_projection_pixel_intensity : Sum of intensity values in projecting
%       pixels in this structure.
% - projection_density : sum_projection_pixels / sum_pixels
% - projection_intensity : sum_projection_pixel_intensity / sum_projection_pixels
% - projection_energy : projection_density * projection_intensity
% - volume : volume of valid pixels in structure. Valid pixels are those
%       not manually annotated as invalid data.
% - projection_volume : volume of projection signal in structure in mm3
% - normalized_projection_volume : projection_volume / total volume of signal
%       in the injection site.
% - max_voxel_density : density of projection signal in 10um3 grid voxel with
%       largest density.
% - max_voxel_x : x (AP, anterior to posterior) coordinate in um of grid voxel with largest density.
% - max_voxel_y : y (DV, dorsal to ventral) coordinate in um of grid voxel with largest density.
% - max_voxel_z : z ("ML", right to left) coordinate in um of grid voxel with largest density.

% Build the URL
url = 'http://connectivity.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json?criteria=[section_data_set_id$eq%d]&num_rows=all';

% Get projection data
status = true;
try
    page = urlread(sprintf(url, experimentID));

catch
    status = false;
    warning('Failed to get data for ID %d\n', experimentID)
end

injectionInfo = struct;
% parse the JSON data
if ~isempty(page)
       
    tmp = jsondecode(page);
    
    if tmp.success
        % load just the injection data 
        injectionInfo = tmp.msg([tmp.msg.is_injection] == 1);
    else
        warning('Failed to jsondecode data for ID %d\n', experimentID)
    end
end

[injectionInfo.experimentID] = deal(experimentID);% store information about the experiment ID

% save results
save([saveFilePath, filesep, 'injectionSummary_all.mat'], 'injectionInfo')

end