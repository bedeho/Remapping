

%
%  AnalyzeDuhamelTruncation.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelTruncation_Result] = AnalyzeDuhamelTruncation(activity, stimuli, stim_control_activity, stim_stimuli)

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
    
    stim_stimuliOffset   = stim_stimuli.stimuliOnset + stim_stimuli.stimuliDuration;
    
    % Analysis process
    responseWindowDuration = 0.200;
    
    % For plots
    saccadeOnsetTimeStep = timeToTimeStep(saccadeOnset, dt);
    responseWindowSize = timeToTimeStep(responseWindowDuration, dt);
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(size(stim_control_activity,3) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');

    %% Latency & Duration
    for p=1:numPeriods,
        
        % Index of R neuron which will recieve remapping activity: currentRF
        truncated_neuronIndex = R_eccentricity + stimuli.stimuli{p}.currentRF + 1;
        
        % Get data for best period of each neuron
        truncate_responseVector  = R_firing_history(truncated_neuronIndex, :, p, 1);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(truncate_responseVector, dt, saccadeOnset, responseWindowDuration);
        
        %% STIM - CONTROL
        
        % Stim response of to offset of stimuli
        stim_responseVector = stim_control_activity(truncated_neuronIndex, :, truncated_neuronIndex);
        stim_stim_offset_response = normalizedIntegration(stim_responseVector, dt, stim_stimuliOffset, responseWindowDuration);
          
        %% FIGURE
        %{
        figure;
       
        remapNumTimeSteps = size(R_firing_history, 2);
        stimNumTimeSteps = size(stim_control_activity, 2);
        
        % Truncate
        subplot(2,2,1);
        hold on;
        rectangle('Position',[saccadeOnsetTimeStep,0,responseWindowSize,1.00],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
        plot(0:(remapNumTimeSteps-1),truncate_responseVector);
        plot([saccadeOnsetTimeStep saccadeOnsetTimeStep],[0 1], 'k');
        title('truncation single neuron response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        xlim([0 (remapNumTimeSteps-1)]);
        
        subplot(2,2,2);
        hold on;
        imagesc(R_firing_history(:, :, p, 1));
        plot([0 (remapNumTimeSteps-1)],[truncated_neuronIndex truncated_neuronIndex], 'w');
        colorbar;
        axis tight;
        colorbar;
        title('remapping population response');
        xlabel(['Time - dt (' num2str(dt) 's']);
            
        % Stim
        subplot(2,2,3);
        hold on;
        offset = timeToTimeStep(stim_stimuliOffset, dt);
        rectangle('Position',[offset,0,responseWindowSize,1.00],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
        plot(0:(stimNumTimeSteps-1),stim_responseVector);
        plot([offset offset],[0 1], 'g');
        title('stimulus control single neuron');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        xlim([0 (stimNumTimeSteps-1)]);
        ylim([0 1]);
        
        subplot(2,2,4);
        hold on;
        imagesc(stim_control_activity(:,:, truncated_neuronIndex));
        plot([1 stimNumTimeSteps],[truncated_neuronIndex truncated_neuronIndex], 'w');
        colorbar;
        axis tight;
        title('stimulus control population response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        %}
        
        %% Save
        
        DuhamelTruncation_Result(p).index                   = truncated_neuronIndex;
        DuhamelTruncation_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        DuhamelTruncation_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        DuhamelTruncation_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelTruncation_Result(p).stim_stim_offset_response = stim_stim_offset_response;
    end
end