function bsv_plotConnectivity3D(injectionSummary, allenAtlasPath, regionToPlot, color, plotPatch)

% load atlas 
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% plot brain outline/volume 
if plotPatch
    [h1, patchHandle] = ya_plotBrainSurface(allenAtlasPath, 'w');
else
    [~, brain_outline] = plotBrainGrid([], []);
end
hold on;

% plot region 
curr_plot_structure_idx = find(contains(st.acronym, regionToPlot));
region_to_plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure_idx(1)}, 2, [])') ./ 255;
slice_spacing = 10;
structure_3d = isosurface(permute(ismember(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end), curr_plot_structure_idx), [3, 1, 2]), 0);
structure_alpha = 0.3;
hold on;
patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', region_to_plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);

% plot injections 
iView = 1;
axis equal

max_voxel_x = injectionSummary.max_voxel_x / 10;
max_voxel_y = injectionSummary.max_voxel_y / 10;
max_voxel_z = injectionSummary.max_voxel_z / 10;

% remove any missing data 
keepMe = max_voxel_x ~= 0;

injection_structure_idx = arrayfun(@(x) find(st.id == injectionSummary.structure_id(x)), 1:size(injectionSummary,1));

plotAllColors = false; %leave false for now
if plotAllColors
    hex_triplets = st.color_hex_triplet(injection_structure_idx(keepMe));
    rgbMatrix = cell2mat(arrayfun(@(x) hex2dec(reshape(x{1}, 2, 3)')', hex_triplets, 'UniformOutput', false)) ./ 255;
else
    rgbMatrix = hex2dec(reshape(st.color_hex_triplet{injection_structure_idx(1)}, 2, [])')' ./ 255;
end
% projection_volume is in mm^ 3, atlas is in 10um x 10um x 10um -> multiply
% by 10^6 - but that looks really odd. QQ leave as 10^3 for now. scatter is
% also the wrong tool for this - it doesn't scale with the plot
dotSize = 10^3;
scatter3(max_voxel_x(keepMe), max_voxel_z(keepMe), max_voxel_y(keepMe), ...
    [injectionSummary.projection_volume(keepMe)].*dotSize , rgbMatrix, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', ...
    1)

set(gcf, 'Color', 'w')
hold off;

% make a movie 
%save as .avi rotating vid
% view([0, 0])
% OptionZ.FrameRate = 5;
% OptionZ.Duration = 15.5;
% OptionZ.Periodic = true;
% CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], 'WellMadeVid', OptionZ)


end
