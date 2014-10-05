
%
%  Sprattling_RemappingPlot_Thesis.m
%  Remapping
%
%  Created by Bedeho Mender on 24/09/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

% Population stim control response plot
function Sprattling_RemappingPlot_Thesis()

    trained_analysis_stim = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune44\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    trained_analysis_remapping = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune44\baseline\TrainedNetwork\analysis-basic-DuhamelRemappingTrace.mat');
    remapping_stimuli = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune44\STIM-basic-DuhamelRemappingTrace\stim.mat');
    
    % Iterate periods
    numPeriods = length(trained_analysis_remapping.DuhamelRemappingTrace_Result);
    RFs = [trained_analysis_stim.Decoded_ReceptiveFieldsLocations];
    
    [B I] = sort(RFs);
    
    num_trials = 17;
    desired_num_cols = 4;
    num_rows = ceil(num_trials / desired_num_cols);
    
    
    period_correct_RF = zeros(1, numPeriods);
    period_predicted_RF = zeros(1, numPeriods);
    
    figure;
    
    for ctr=1:numPeriods,
        
        subplot(num_rows, desired_num_cols, ctr);
        
        correct_RF = remapping_stimuli.stimuli{ctr}.headCenteredTargetLocations-remapping_stimuli.stimuli{ctr}.saccadeTargets;
        period_correct_RF(ctr) = correct_RF;
        
        RIs = trained_analysis_remapping.DuhamelRemappingTrace_Result(ctr).remapping_index_all;
        
        predicted_RF = dot(RFs, RIs)/sum(RIs);
        period_predicted_RF(ctr) = predicted_RF;
        
        
        plot(RFs, RIs, 'o', 'MarkerSize', 4);
        hold on;
        plot([correct_RF correct_RF],[0  sqrt(2)],'-g');
        hold on;
        plot([predicted_RF predicted_RF],[0  sqrt(2)],'-b');
 

        ylim([0 sqrt(2)]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
    
    p = polyfit(period_correct_RF, period_predicted_RF,1)
    lin = p(1)*period_correct_RF + p(2);
    lin_straight = period_correct_RF + p(2);
    
    
    figure;
    plot(period_correct_RF, period_predicted_RF , 'o');
    hold on;
    plot(period_correct_RF, lin ,'b');
    plot(period_correct_RF, lin_straight ,'g');
    xlabel('Correct Retinal Location (deg)');
    ylabel('Decoded Retinal Location (deg)');
    xlim([-45 45]);
    ylim([-45 45]);
    corrcoef(period_correct_RF,period_predicted_RF)
    
    
end