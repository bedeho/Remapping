
%
%  AnalyzeLHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [LHeiserAnalysis] = AnalyzeLHeiser(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli)

    % Do underlying analysis
    LHeiser_DuhamelRemappingTrace = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli);
    
    % Save indexes for all neurons and directions in 2d-array
    Unique_Training_RF_Locations = unique(stimuli.currentRF);
    numUnique = length(Unique_Training_RF_Locations);
    numberOfDirections = activity.numPeriods/numUnique;
    uniqueIndexes = reshape([LHeiser_DuhamelRemappingTrace.remapping_index], numberOfDirections, numUnique);
    
    % Count num directions per neuron
    uniqueResponseCount = sum(uniqueIndexes > 0.1);
    
    %% Selectivity
     
    
    
    % Save
    LHeiserAnalysis.numUnique                    = numUnique;
    LHeiserAnalysis.numberOfDirections           = numberOfDirections;
    LHeiserAnalysis.uniqueIndexes                = uniqueIndexes;
    LHeiserAnalysis.Unique_Training_RF_Locations = Unique_Training_RF_Locations;
    LHeiserAnalysis.uniqueResponseCount          = uniqueResponseCount;
    
end