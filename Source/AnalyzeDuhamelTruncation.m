
%
%  AnalyzeDuhamelTruncation.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelTruncation_Result] = AnalyzeDuhamelTruncation(activity, stimuli, StimuliControl_Result)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        %LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-DuhamelTruncation.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-DuhamelTruncation/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    saccadeOnset         = stimuli.saccadeOnset;
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(length(StimuliControl_Result) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');

    
    %% Latency & Duration
    for p=1:numPeriods,
        
        % Index of R neuron which will recieve remapping activity: currentRF
        remappedInto_neuronIndex = R_eccentricity + stimuli.stimuli{p}.currentRF + 1;
        
        % Get data for best period of each neuron
        remap_responseVector  = R_firing_history(remappedInto_neuronIndex, :, p, 1);

        % Find latency and duration
        [latencyTimeStep, duration] = findNeuronalLatency(responseThreshold, remap_responseVector, latencyWindowLength);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(remap_responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(remap_responseVector, dt, stimuliOnset, responseWindowDuration);
        
        %% STIM - CONTROL
        
        % Index of R neuron which has RF over stim location before saccade: futureRF
        futureRF_neuronIndex = R_eccentricity + stimuli.stimuli{p}.stim_screen_location + 1;
        
        % Stim response of remapped neuron when stim is in future RF
        stim_responseVector = stim_control_activity(remappedInto_neuronIndex, :, futureRF_neuronIndex);
        stim_control_response = normalizedIntegration(stim_responseVector, dt, stim_stimuliOnset + stim_responseWindowStart, responseWindowDuration);
        stim_index = saccadeonset_response - stim_control_response;
        
        % Find latency and duration
        stim2_responseVector = stim_control_activity(futureRF_neuronIndex, :, futureRF_neuronIndex);
        [stim_latencyTimeStep stim_duration] = findNeuronalLatency(responseThreshold, stim2_responseVector, latencyWindowLength);
        
        %% SACC = CONTROL
        
        % Index of R neuron which has RF representing the saccade
        % performed
        sacc_neuronIndex = S_eccentricity + stimuli.stimuli{p}.saccadeTargets + 1;
        
        % Sacc response of remapped neuron when given saccade is performed
        sacc_responseVector = sacc_control_activity(remappedInto_neuronIndex, :, sacc_neuronIndex);
        saccade_control_response = normalizedIntegration(sacc_responseVector, dt, sacc_saccadeOnset, responseWindowDuration);
        sacc_index = saccadeonset_response - saccade_control_response;
        
        % Remapping Indexes
        remapping_index = sqrt(stim_index^2 + sacc_index^2); 
        
        %% FIGURE
        %{
        figure;
       
        remapNumTimeSteps = size(R_firing_history, 2);
        stimNumTimeSteps = size(stim_control_activity, 2);
        saccNumTimeSteps = size(sacc_control_activity,2);
        
        % Remap
        subplot(3,2,1);
        hold on;
        rectangle('Position',[saccadeOnsetTimeStep,0,responseWindowSize,1.00],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
        plot(0:(remapNumTimeSteps-1),remap_responseVector);
        plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        plot([stimuliOnsetTimeStep stimuliOnsetTimeStep],[0 1], 'g');
        plot([saccadeOnsetTimeStep saccadeOnsetTimeStep],[0 1], 'k');
        title('remapping single neuron response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        xlim([0 (remapNumTimeSteps-1)]);
        
        subplot(3,2,2);
        hold on;
        imagesc(R_firing_history(:, :, p, 1));
        plot([0 (remapNumTimeSteps-1)],[remappedInto_neuronIndex remappedInto_neuronIndex], 'w');
        colorbar;
        axis tight;
        colorbar;
        title('remapping population response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        
        %}
        
        %% Save
        %{
        DuhamelTruncation_Result(p).index                   = remappedInto_neuronIndex;
        DuhamelTruncation_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        DuhamelTruncation_Result(p).futureRF                = stimuli.stimuli{p}.stim_screen_location;
        DuhamelTruncation_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        DuhamelTruncation_Result(p).remappingLatency        = stepToTime(latencyTimeStep, dt)-saccadeOnset;
        DuhamelTruncation_Result(p).stimLatency             = stepToTime(stim_latencyTimeStep, dt)-stim_stimuliOnset;
        
        DuhamelTruncation_Result(p).Duration                = duration*dt;
        DuhamelTruncation_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelTruncation_Result(p).stimulionset_response   = stimulionset_response;
        
        DuhamelTruncation_Result(p).stim_index              = stim_index;
        DuhamelTruncation_Result(p).sacc_index              = sacc_index;
        DuhamelTruncation_Result(p).remapping_index         = remapping_index;
        %}
        
    end
   
    
    
end