
%
%  Sprattling_StimPlot_Thesis.m
%  Remapping
%
%  Created by Bedeho Mender on 24/09/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

% Population stim control response plot
function Sprattling_StimPlot_Thesis()

    untrained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune52\baseline\BlankNetwork\analysis-basic-StimuliControl.mat');
    trained_analysis = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune52\baseline\TrainedNetwork\analysis-basic-StimuliControl.mat');
    trained_analysis_remapping = load('C:\Users\bedeho\Documents\GitHub\Remapping\Experiments\sprattling_visual_learning_bigepoch20-connectivitfix-tune52\baseline\TrainedNetwork\analysis-basic-DuhamelRemappingTrace.mat');
 
    choosen_neurons = [trained_analysis_remapping.DuhamelRemappingTrace_Result.index];
    
    num_neurons = 91; %stim_control_activity.R_N;
    desired_num_cols = 10;
    num_rows = ceil(num_neurons / desired_num_cols);

    % Iterate neurons
    for ctr=1:num_neurons,
        
        subplot(num_rows, desired_num_cols, ctr);

        plot(untrained_analysis.cross_trial_activity(ctr, :) ,'r');
        hold on;
        plot(trained_analysis.cross_trial_activity(ctr, :) ,'b');
        hold on;
        plot(get_x(trained_analysis.Decoded_ReceptiveFieldsLocations(ctr)), [0 1],'b');
        
        %{
        number_of_times_this_neuron_was_picked = nnz(choosen_neurons == ctr);
        
        if(number_of_times_this_neuron_was_picked == 1),
            text(1,0.8,num2str(find(choosen_neurons == ctr)),'FontSize',6)
        elseif(number_of_times_this_neuron_was_picked > 1),
            warning('one neuron picked more than ones, is that ok??');
        end 
        %}
        
        ylim([0 1]);
        xlim([1 num_neurons]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
    
    function x = get_x(loc)
        x = (loc + 46) * [1 1];
    end
    
    
end