function nsv_plotConnectivity3D(injectionSummary, allenAtlasPath, region, color)

tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean


if plotPatch
    ya_plotBrainSurface(allenAtlasPath)
else
    [~, brain_outline] = plotBrainGrid([], []);
end
hold on;

% injection coordinates
combinedStr = strcat(injectionTable.injection_coordinates{:});
combinedStr = strrep(combinedStr, ']', ', ');
cleanStr = erase(combinedStr, {'[', ']'});
cleanStr(end-1:end) = '';
numStrs = strsplit(cleanStr, ',');
nums = str2double(numStrs);
injection_coordinatesMatrix = reshape(nums, 3, []).'; % Adjust dimensions as necessary

% Convert cell array to char array for vectorized operations
hexMatrix = char(injectionTable .structure_color);
rHex = hexMatrix(:, 1:2);
gHex = hexMatrix(:, 3:4);
bHex = hexMatrix(:, 5:6);
r = hex2dec(rHex);
g = hex2dec(gHex);
b = hex2dec(bHex);
rgbMatrix = [r, g, b];


theseInjections = ismember(injectionTable.structure_abbrev, injectionAreas(1:end-1));

curr_plot_structure_idx = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = regionColors{1};

structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure_idx, [3, 1, 2]), 0);
structure_alpha = 0.2;


[h1, patchHandle] = ya_plotBrainSurface(allenAtlasPath, 'w');


hold on;
structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', plot_structure_color, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);

iView = 1;
axis equal
[ap_max, dv_max, ml_max] = size(tv);

hold on;


max_voxel_x = [injection_coordinatesMatrix(theseInjections, 1)] / 10;
max_voxel_y = [injection_coordinatesMatrix(theseInjections, 2)] / 10;
max_voxel_z = [injection_coordinatesMatrix(theseInjections, 3)] / 10;

structure_color = rgbMatrix(theseInjections, :) ./ 255;

scatter3(max_voxel_x, max_voxel_z, max_voxel_y, ...
    [sqrt(injectionTable.injection_volume(theseInjections))]*50, structure_color, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', ...
    1)
if iView ~= 1
    set(h1, 'Ydir', 'reverse')
end

set(gcf, 'Color', 'w')


end
