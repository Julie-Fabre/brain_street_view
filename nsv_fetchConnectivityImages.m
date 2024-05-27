function status = nsv_fetchConnectivityImages(experimentID, saveFilePath)
% density.raw is in AP x DV x ML, 100um resolution
status = true;
zipFile = [saveFilePath, filesep, 'temp.zip'];
urlwrite(sprintf('http://api.brain-map.org/grid_data/download/%d?include=density', experimentID), zipFile);
unzip(zipFile, num2str(experimentID));

filePath = [saveFilePath, '/density.raw'];
if ~exist(filePath, 'file')
    status=false;
    warning('Error downloading experiment %s, skipping.', experimentID);
end

end