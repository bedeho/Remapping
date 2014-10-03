
%
%  Sprattling_RemappingPlot_Thesis.m
%  Remapping
%
%  Created by Bedeho Mender on 24/09/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

% Population stim control response plot
function Sprattling_RemappingPlot_Thesis()

    trained_analysis_stim = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune19\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    trained_analysis_remapping = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune19\baseline\TrainedNetwork\analysis-basic-DuhamelRemappingTrace.mat');
   
    %trained_analysis_stim = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\classic\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    %trained_analysis_remapping = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\classic\baseline\TrainedNetwork\analysis-basic-DuhamelRemappingTrace.mat');
   
    remapping_stimuli = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix\STIM-basic-DuhamelRemappingTrace\stim.mat');
    
    % Iterate periods
    numPeriods = length(trained_analysis_remapping.DuhamelRemappingTrace_Result);
    RFs = [trained_analysis_stim.Decoded_ReceptiveFieldsLocations];
    
    [B I] = sort(RFs);
    
    for ctr=1:numPeriods,
        
        figure;
        
        correct_RF = remapping_stimuli.stimuli{ctr}.headCenteredTargetLocations-remapping_stimuli.stimuli{ctr}.saccadeTargets
        
        % FIX LATER
        %RIs = trained_analysis_remapping.DuhamelRemappingTrace_Result(ctr).remapping_index_all;
        % FIX LATER
        
        RIs = trained_analysis_remapping.DuhamelRemappingTrace_Result(ctr).remapping_response_all; 
        
        plot(RFs, RIs, 'o');
        hold on;
        plot([correct_RF correct_RF],[0  sqrt(2)],'-r');
 

        ylim([0 sqrt(2)]);
    end
    
end