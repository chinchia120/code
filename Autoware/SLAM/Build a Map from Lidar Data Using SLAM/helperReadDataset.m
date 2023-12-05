%% =============== Function of helperReadDataset =============== %%
function datasetTable = helperReadDataset(dataFolder, pointCloudFilePattern)
    % Create a file datastore to read in files in the right order
    fileDS = fileDatastore(pointCloudFilePattern, 'ReadFcn', @helperReadPointCloudFromFile);

    % Extract the file list from the datastore
    pointCloudFiles = fileDS.Files;

    imuConfigFile = fullfile(dataFolder, 'scenario1', 'imu.cfg');
    insDataTable = helperReadINSConfigFile(imuConfigFile);

    % Delete the bad row from the INS config file
    insDataTable(1447, :) = [];

    % Remove columns that will not be used
    datasetTable = removevars(insDataTable, ...
        {'Num_Satellites', 'Latitude', 'Longitude', 'Altitude', 'Omega_Heading', ...
         'Omega_Pitch', 'Omega_Roll', 'V_X', 'V_Y', 'V_ZDown'});

    datasetTable = addvars(datasetTable, pointCloudFiles, 'Before', 1, 'NewVariableNames', "PointCloudFileName");
end