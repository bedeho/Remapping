
%
%  SaveActivity.m
%  Remapping
%
%  Created by Bedeho Mender on 05/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function SaveActivity( fileName , R_N ...
                                , S_N ...
                                , C_N ...
                                , numEpochs ...
                                , numPeriods ...
                                , dt ...
                                , V_firing_history ...
                                , R_firing_history ...
                                , S_firing_history ...
                                , C_firing_history ...
                                , V_activation_history ...
                                , R_activation_history ...
                                , S_activation_history ...
                                , C_activation_history)
                            
    % Open file
    fileID = fopen(fileName, 'w');
    
    % Write to file
    fwrite(fileID, R_N, 'uint32');
    fwrite(fileID, S_N, 'uint32');
    fwrite(fileID, C_N, 'uint32');
    fwrite(fileID, numEpochs, 'uint32');
    fwrite(fileID, numPeriods, 'uint32');
    
    fwrite(fileID, dt, 'float32');
    fwrite(fileID, V_firing_history, 'float32');
    fwrite(fileID, R_firing_history, 'float32');
    fwrite(fileID, S_firing_history, 'float32');
    fwrite(fileID, C_firing_history, 'float32');
    fwrite(fileID, V_activation_history, 'float32');
    fwrite(fileID, R_activation_history, 'float32');
    fwrite(fileID, S_activation_history, 'float32');
    fwrite(fileID, C_activation_history, 'float32');
    
    % Close file
    fclose(fileID);
    
end