
%
%  AnalyzeLHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [LHeiserAnalysis] = AnalyzeLHeiser(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli)

    % Do underlying analysis
    LHeiser_DuhamelRemappingTrace = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli);

    numPeriods           = activity.numPeriods;
    

    for p=1:numPeriods,
        
        stimuli.stimuli{p}.currentRF
        
        stimuli.stimuli{p}.futureRF
        
        stimuli.stimuli{p}.saccadeTargets
        
  
    end
    
end