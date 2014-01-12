
%
%  AnalyzeDuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelRemappin_Result] = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        %LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-DuhamelRemapping.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-DuhamelRemapping/stim.mat');
    end
    
    % Get data remapping activity data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    S_eccentricity       = stimuli.S_eccentricity;
    saccadeOnset         = stimuli.saccadeOnset;
    stimuliOnset         = stimuli.stimuliOnset;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    
    stim_stimuliOnset = stim_stimuli.stimuliOnset;
    sacc_saccadeOnset = sacc_stimuli.saccadeOnset;
    
    % Analysis params
    %latencyWindowSize   = 0.020; % (s), colby papers
    %latencyWindowLength = timeToTimeStep(latencyWindowSize, dt);
    responseWindowDuration = 0.200;
    
    stim_responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    
    % For plots
    stimuliOnsetTimeStep = timeToTimeStep(stimuliOnset, dt);
    saccadeOnsetTimeStep = timeToTimeStep(saccadeOnset, dt);
    stim_stimuliOnsetTimeStep = timeToTimeStep(stim_stimuliOnset, dt);
    sacc_saccadeOnsetTimeStep = timeToTimeStep(sacc_saccadeOnset, dt);
    
    responseWindowSize = timeToTimeStep(responseWindowDuration, dt);
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(size(stim_control_activity,3) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');
    assert(size(sacc_control_activity,3) == 2*S_eccentricity+1, 'Saccade control must be done for saccade map');
    
    %% Latency & Duration
    for p=1:numPeriods,
        
        % Index of R neuron which will recieve remapping activity: currentRF
        remappedInto_neuronIndex = R_eccentricity + stimuli.stimuli{p}.currentRF + 1;
        
        % Get data for neuron
        remap_responseVector  = R_firing_history(remappedInto_neuronIndex, :, p, 1);

        % Find latency and duration
        %[latencyTimeStep, duration] = findNeuronalLatency(remap_responseVector, latencyWindowLength);
        latencyTimeStep = findNeuronalLatency_NEW(remap_responseVector, dt);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(remap_responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(remap_responseVector, dt, stimuliOnset, responseWindowDuration);
        
        %% STIM - CONTROL
        
        % Index of R neuron which has RF over stim location before saccade: futureRF
        futureRF_neuronIndex = R_eccentricity + stimuli.stimuli{p}.futureRF + 1;
        
        % Stim response of remapped neuron when stim is in future RF
        stim_responseVector = stim_control_activity(remappedInto_neuronIndex, :, futureRF_neuronIndex);
        stim_control_response = normalizedIntegration(stim_responseVector, dt, stim_stimuliOnset + stim_responseWindowStart, responseWindowDuration);
        stim_index = saccadeonset_response - stim_control_response;
        
        % Find latency and duration
        stim2_responseVector = stim_control_activity(futureRF_neuronIndex, :, futureRF_neuronIndex);
        %[stim_latencyTimeStep stim_duration] = findNeuronalLatency(stim2_responseVector, latencyWindowLength);
        stim_latencyTimeStep= findNeuronalLatency_NEW(stim2_responseVector, dt);
        
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
       
        R_N = size(R_firing_history, 1);
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
        plot([saccadeOnsetTimeStep saccadeOnsetTimeStep],[0 1], '--k');
        title('remapping single neuron response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        xlim([0 (remapNumTimeSteps-1)]);
        
        subplot(3,2,2);
        hold on;
        imagesc(R_firing_history(:, :, p, 1));
        plot([0 (remapNumTimeSteps-1)],[remappedInto_neuronIndex remappedInto_neuronIndex], 'w');
        plot([latencyTimeStep latencyTimeStep],[1 R_N], 'r');
        plot([stimuliOnsetTimeStep stimuliOnsetTimeStep],[1 R_N], 'g');
        plot([saccadeOnsetTimeStep saccadeOnsetTimeStep],[1 R_N], '--k');
        colorbar;
        axis tight;
        colorbar;
        title('remapping population response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        
        % Stim
        subplot(3,2,3);
        hold on;
        onset = timeToTimeStep(stim_stimuliOnset + stim_responseWindowStart, dt);
        rectangle('Position',[onset,0,responseWindowSize,1.00],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
        plot(0:(stimNumTimeSteps-1),stim_responseVector);
        plot([stim_stimuliOnsetTimeStep stim_stimuliOnsetTimeStep],[0 1], 'g');
        title('stimulus control single neuron');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        xlim([0 (stimNumTimeSteps-1)]);
        ylim([0 1]);
        
        subplot(3,2,4);
        hold on;
        imagesc(stim_control_activity(:,:, futureRF_neuronIndex));
        plot([1 stimNumTimeSteps],[remappedInto_neuronIndex remappedInto_neuronIndex], 'w');
        colorbar;
        axis tight;
        title('stimulus control population response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        
        % Sacc
        subplot(3,2,5);
        hold on;
        
        rectangle('Position',[sacc_saccadeOnsetTimeStep,0,responseWindowSize,1.00],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
        
        plot(0:(saccNumTimeSteps-1), sacc_responseVector);
        plot([sacc_saccadeOnsetTimeStep sacc_saccadeOnsetTimeStep],[0 1], 'g');
        title('saccade control single neuron');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        ylim([0 1]);
        
        subplot(3,2,6);
        hold on;
        imagesc(sacc_control_activity(:, :, sacc_neuronIndex));
        plot([1 saccNumTimeSteps],[remappedInto_neuronIndex remappedInto_neuronIndex], 'w');
        colorbar;
        axis tight;
        title('saccade population response');
        xlabel(['Time - dt (' num2str(dt) ' s)']);
        

        %}
        
        %% Save
        DuhamelRemappin_Result(p).index                   = remappedInto_neuronIndex;
        DuhamelRemappin_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        DuhamelRemappin_Result(p).futureRF                = stimuli.stimuli{p}.futureRF;
        DuhamelRemappin_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        DuhamelRemappin_Result(p).remappingLatency        = (latencyTimeStep - timeToTimeStep(saccadeOnset,dt))*dt;
        DuhamelRemappin_Result(p).stimLatency             = (stim_latencyTimeStep - timeToTimeStep(stim_stimuliOnset, dt))*dt;
        
        %DuhamelRemappin_Result(p).Duration                = duration*dt;
        DuhamelRemappin_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelRemappin_Result(p).stimulionset_response   = stimulionset_response;
        
        DuhamelRemappin_Result(p).stim_index              = stim_index;
        DuhamelRemappin_Result(p).sacc_index              = sacc_index;
        DuhamelRemappin_Result(p).remapping_index         = remapping_index;
        
    end
    
end