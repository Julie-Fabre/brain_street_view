function nsv_fetchConnectivityImages(experimentID, saveFilePath)
% density.raw is in AP x DV x ML, 100um resolution

zipFile = [saveFilePath, filesep, 'temp.zip'];
urlwrite(sprintf('http://api.brain-map.org/grid_data/download/%d?include=density', experimentID), zipFile);
unzip(zipFile, num2str(experimentID));

filePath = [saveFilePath, '/density.raw'];
if ~exist(filePath, 'file')
    error('Error downloading experiment %s .', experimentID);
end

end