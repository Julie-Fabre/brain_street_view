function status = fetchConnectivityImages(experimentID, saveFilePath)
% density.raw is in AP x DV x ML, 100um resolution
status = true;

% Create directory if it doesn't exist
if ~exist(saveFilePath, 'dir')
    [mkdirStatus, msg] = mkdir(saveFilePath);
    if ~mkdirStatus
        status = false;
        warning('Failed to create directory %s: %s', saveFilePath, msg);
        return;
    end
end

% Verify directory exists and is writable
if ~exist(saveFilePath, 'dir')
    status = false;
    warning('Directory %s does not exist after creation attempt', saveFilePath);
    return;
end

% Additional check - try to create a test file to verify write permissions
testFile = fullfile(saveFilePath, 'test_write.txt');
try
    fid = fopen(testFile, 'w');
    if fid == -1
        status = false;
        warning('Cannot write to directory %s - check permissions', saveFilePath);
        return;
    end
    fclose(fid);
    delete(testFile);
catch
    status = false;
    warning('Cannot write to directory %s - check permissions', saveFilePath);
    return;
end

zipFile = [saveFilePath, filesep, 'temp.zip'];

% Try to download the file with error handling
try
    % First try urlwrite
    urlwrite(sprintf('http://api.brain-map.org/grid_data/download/%d?include=density', experimentID), zipFile);
catch ME
    % If urlwrite fails, try websave as a fallback (newer MATLAB versions)
    try
        websave(zipFile, sprintf('http://api.brain-map.org/grid_data/download/%d?include=density', experimentID));
    catch ME2
        status = false;
        warning('Failed to download data for experiment %d: %s', experimentID, ME.message);
        return;
    end
end

% Unzip the file
try
    unzip(zipFile, saveFilePath);
catch ME
    status = false;
    warning('Failed to unzip data for experiment %d: %s', experimentID, ME.message);
    return;
end

filePath = [saveFilePath, '/density.raw'];
if ~exist(filePath, 'file')
    status = false;
    warning('Error downloading experiment %s, skipping.', experimentID);
end

end