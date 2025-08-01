
% get striatum subdivisions described in Hintyrian / Dong Nat Neuro 2017
% implement in Allen Atlas by Yongsoo Kim (Chon et al., 2019, Nature
% Communications (PMID: 31699990), DOI: 10.1038/s41467-019-13057-w)
%   - integrated detailed CP segmentation into the Allen CCF.
% The latest file can be downloaded at https://figshare.com/articles/dataset/Unified_mouse_brain_atlas_v2/25750983
% - Here is the related X post from my lab account. 
% https://twitter.com/yongsookimlab/status/1786782412397011372


addpath('/path/to/NIfTI_20140122');  % Add path to NIfTI toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
addpath('/path/to/npy-matlab');      % Add path to npy-matlab toolbox

% Specify input and output file paths
input_file = '/home/julie/Dropbox/Atlas/allenCCF_v2/UnifiedAtlas_Label_v2_20um-isotropic.nii';
output_file = '/home/julie/Dropbox/Atlas/allenCCF_v2/annotation_volume_v2_20um_by_index.npy';

% Load NIfTI file
nii = load_nii(input_file);

% Extract data from NIfTI structure
data = nii.img;

% Re-order from DV x ML x AP to AP x ML x DV
data_pa_ml_dv_20um = permute(data, [3, 2, 1]);

% reorder from posterior-to-anterior to anterior-to-poster 
data_ap_ml_dv_20um = data_pa_ml_dv_20um(660:-1:1,:,:);

figure(); imagesc(squeeze(data_ap_ml_dv_20um((810/2),:,:))) % visualize 

% Save data as .npy file
writeNPY(data_ap_ml_dv_20um, output_file); 

%% GPe and SNr, recreate labels from Foster/Dond Nature 2021 
foster_path = '/home/julie/Dropbox/MATLAB/onPaths/basalganglia/static/files/';
%load in a slice, use indexes to replace. if none: nearest value. 
snr_81_table = readtable([foster_path, 'SNr81_boxgrid_data.csv']);


allenv1_81=imread('/home/julie/Downloads/100142143_322 (1).jpg');
size(allenv1)
%% prev gpe snr tries
% Original matrix dimensions and pixel size
original_dims = [660, 400, 570];
original_pixel_size = 20; % in micrometers

% New pixel size
new_pixel_size = 63; % in micrometers

% Calculate scaling factor
scale_factor = original_pixel_size / new_pixel_size;

% Calculate new dimensions (rounding down)
new_dims = [660, floor(original_dims(2:3) * scale_factor)];

% Create example data (replace this with your actual data)
original_matrix = data_ap_ml_dv_20um;

% Perform the transformation
transformed_matrix = imresize3(data_ap_ml_dv_20um, new_dims, 'linear');
smol = squeeze(transformed_matrix(810/2+1, 75+3:90+3, 55+3:66+6));
%% ara 81
load('snr_ara81.mat')
snr_foster_81(snr_foster_81>0) = snr_foster_81(snr_foster_81>0) + 2000;

transformed_matrix_snr_81_init = squeeze(transformed_matrix(810/2,:,:));
transformed_matrix_snr_81 = transformed_matrix_snr_81_init;
transformed_matrix_snr_81(transformed_matrix_snr_81_init~=snr_id)=0;

transformed_matrix_snr_81(transformed_matrix_snr_81_init~=snr_id) = transformed_matrix_snr_81_init(transformed_matrix_snr_81_init~=snr_id);

transformed_matrix(810/2+1,:,:) = transformed_matrix_snr_81;

%% ara 83
%snr_foster_83 = 
%%load('snr_ara81.mat')
%snr_foster_81(snr_foster_81>0) = snr_foster_81(snr_foster_81>0) + 2000;

transformed_matrix_snr_91_init = squeeze(transformed_matrix(910/2,:,:));
transformed_matrix_snr_91 = transformed_matrix_snr_91_init;
transformed_matrix_snr_91(transformed_matrix_snr_91_init~=snr_id)=0;
transformed_matrix_91 = transformed_matrix_snr_91;
figure();imagesc(transformed_matrix_snr_91)
clim([2000, 2020])
transformed_matrix_snr_91(transformed_matrix_snr_91_init~=snr_id) = transformed_matrix_snr_91_init(transformed_matrix_snr_91_init~=snr_id);

%transformed_matrix(830/2+1,:,:) = transformed_matrix_snr_83;
save('full_ara91.mat', 'transformed_matrix_snr_91');