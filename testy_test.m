cl_myPaths;

plotType = '2D';
plotInjections = true;
nChunks = 10;
targetColors = rgb('HotPink');
injectionColors = rgb('Gold');
test_expID = 590987294;
nsv_getAllenAtlasProjections(test_expID, {'CP'}, allenAtlasPath, localDataPath, plotType, plotInjections, nChunks, targetColors, injectionColors)

%allen atlas
cl_getAllenAtlasProjections_visual2('VIS_test', '2D', 10)


%allen atlas
cl_getAllenAtlasProjections_visual3('VIS_test', '2D', 10)
%allen atlas
cl_getAllenAtlasProjections_visual3('VIS', '2D', 10)

cl_getAllenAtlasProjections_visual_clean('VIS_test', '2D', 10)
cl_getAllenAtlasProjections_visual_clean('VIS', '2D', 10)
cl_getAllenAtlasProjections_visual_clean_imgs('VIS', '2D', 10)
%brainreg 
cl_getAllenAtlasProjections_visual('VIS_test', '2D', 10)

