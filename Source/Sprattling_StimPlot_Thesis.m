
%
%  Sprattling_StimPlot_Thesis.m
%  Remapping
%
%  Created by Bedeho Mender on 24/09/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

% Population stim control response plot
function Sprattling_StimPlot_Thesis()

    declareGlobalVars();
    global THESIS_FIGURE_PATH;

    %untrained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20\baseline\BlankNetwork\analysis-basic-StimuliControl.mat');
    %trained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    
    untrained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix\baseline\BlankNetwork\analysis-basic-StimuliControl.mat');
    trained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    
    stimuli = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix\STIM-basic-StimuliControl\stim.mat');
    
    num_neurons = 91; %stim_control_activity.R_N;
    desired_num_cols = 10;
    num_rows = ceil(num_neurons / desired_num_cols);

    % Iterate neurons
    for ctr=1:num_neurons,
        
        subplot(num_rows, desired_num_cols, ctr);

        plot(untrained_analysis.cross_trial_activity(ctr, :) ,'r');
        hold on;
        plot(trained_analysis.cross_trial_activity(ctr, :) ,'b');

        ylim([0 1]);
        xlim([1 num_neurons]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
    
    
end