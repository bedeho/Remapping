
%
%  LoadActivity.m
%  Remapping
%
%  Created by Bedeho Mender on 05/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function activity = LoadActivity(fileName)
                            
    % Open file
    fileID = fopen(fileName, 'w');
    
    % Write to file
    activity.R_N                  = fread(fileID, 'uint32');
    activity.S_N                  = fread(fileID, 'uint32');
    activity.C_N                  = fread(fileID, 'uint32');
    activity.numEpochs            = fread(fileID, 'uint32');
    activity.numPeriods           = fread(fileID, 'uint32');
    
    activity.dt                   = fread(fileID, 'float32');
    
    activity.V_firing_history     = fread(fileID, [activity.R_N activity.numPeriods], 'float32');
    activity.R_firing_history     = fread(fileID, [activity.R_N activity.numPeriods], 'float32');
    activity.S_firing_history     = fread(fileID, [activity.S_N activity.numPeriods], 'float32');
    activity.C_firing_history     = fread(fileID, [activity.C_N activity.numPeriods], 'float32');
    
    activity.V_activation_history = fread(fileID, [activity.R_N activity.numPeriods], 'float32');
    activity.R_activation_history = fread(fileID, [activity.R_N activity.numPeriods], 'float32');
    activity.S_activation_history = fread(fileID, [activity.S_N activity.numPeriods], 'float32');
    activity.C_activation_history = fread(fileID, [activity.C_N activity.numPeriods], 'float32');
    
    % Close file
    fclose(fileID);
    
end