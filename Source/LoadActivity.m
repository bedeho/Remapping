
%
%  LoadActivity.m
%  Remapping
%
%  Created by Bedeho Mender on 05/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function activity = LoadActivity(fileName)
                            
    % Open file
    fileID = fopen(fileName, 'r');
    
    % Write to file
    activity.R_N                  = fread(fileID, 1, 'uint32');
    activity.S_N                  = fread(fileID, 1, 'uint32');
    activity.C_N                  = fread(fileID, 1, 'uint32');
    activity.numEpochs            = fread(fileID, 1, 'uint32');
    activity.numPeriods           = fread(fileID, 1, 'uint32');
    activity.numSavedTimeSteps    = fread(fileID, 1, 'uint32');
    activity.outputSavingRate     = fread(fileID, 1, 'uint32');
    
    activity.dt                   = fread(fileID, 1, 'float32');
    

    size = [activity.numSavedTimeSteps activity.numPeriods activity.numEpochs];

    
    activity.V_firing_history     = fread(fileID, [activity.R_N size], 'float32');
    activity.R_firing_history     = fread(fileID, [activity.R_N size], 'float32');
    activity.S_firing_history     = fread(fileID, [activity.S_N size], 'float32');
    activity.C_firing_history     = fread(fileID, [activity.C_N size], 'float32');
    
    activity.V_activation_history = fread(fileID, [activity.R_N size], 'float32');
    activity.R_activation_history = fread(fileID, [activity.R_N size], 'float32');
    activity.S_activation_history = fread(fileID, [activity.S_N size], 'float32');
    activity.C_activation_history = fread(fileID, [activity.C_N size], 'float32');
    
    % Close file
    fclose(fileID);
    
end