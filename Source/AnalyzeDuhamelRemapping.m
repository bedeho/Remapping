
%
%  AnalyzeDuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelRemappin_Result] = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity, sacc_control_activity)

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
    
    % Analysis params
    latencyWindowSize   = 0.020; % (s), colby papers
    latencyWindowLength = timeToTimeStep(latencyWindowSize, dt);
    responseWindowDuration = 0.200;
    responseThreshold   = 0.5;
    
    responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(size(stim_control_activity,3) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');
    assert(size(sacc_control_activity,3) == 2*S_eccentricity+1, 'Saccade control must be done for saccade map');
    
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
        
        %% Control activity
        
        % Index of R neuron which has RF over stim location before saccade: futureRF
        futureRF_neuronIndex = R_eccentricity + stimuli.stimuli{p}.stim_screen_location + 1;
        
        % Stim response of remapped neuron when stim is in future RF
        stim_responseVector = stim_control_activity(remappedInto_neuronIndex, :, futureRF_neuronIndex);
        stim_control_response = normalizedIntegration(stim_responseVector, dt, stimuliOnset + responseWindowStart, responseWindowDuration);
        stim_index = saccadeonset_response - stim_control_response;
        
        % Index of R neuron which has RF representing the saccade
        % performed
        sacc_neuronIndex = R_eccentricity + stimuli.stimuli{p}.saccadeTargets + 1;
        
        % Sacc response of remapped neuron when given saccade is performed
        sacc_responseVector = sacc_control_activity(remappedInto_neuronIndex, :, sacc_neuronIndex);
        saccade_control_response = normalizedIntegration(sacc_responseVector, dt, saccadeOnset, responseWindowDuration);
        sacc_index = saccadeonset_response - saccade_control_response;
        
        % Remapping Indexes
        remapping_index = sqrt(stim_index^2 + sacc_index^2); 
        
        %% CONTROL FIGURE
        figure;
        
        % Remap
        subplot(3,2,1); hold on;
        plot(remap_responseVector);
        plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        title('remapping single neuron response');
        
        subplot(3,2,2);
        imagesc(R_firing_history(:, :, p, 1)); axis tight;
        title('remapping population response');
        
        % Stim
        subplot(3,2,3);
        plot(stim_responseVector);
        subplot(3,2,4);
        imagesc(stim_control_activity(:,:, futureRF_neuronIndex)); axis tight;
        
        % Sacc
        subplot(3,2,5);
        plot(sacc_responseVector);
        subplot(3,2,6);
        imagesc(sacc_control_activity(:, :, sacc_neuronIndex)); axis tight;
        
        %% Save
        DuhamelRemappin_Result(p).index                   = remappedInto_neuronIndex;
        DuhamelRemappin_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        DuhamelRemappin_Result(p).futureRF                = stimuli.stimuli{p}.stim_screen_location;
        DuhamelRemappin_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        
        DuhamelRemappin_Result(p).latency                 = stepToTime(latencyTimeStep, dt)-saccadeOnset;
        DuhamelRemappin_Result(p).Duration                = duration*dt;
        DuhamelRemappin_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelRemappin_Result(p).stimulionset_response   = stimulionset_response;
        
        DuhamelRemappin_Result(p).stim_index              = stim_index;
        DuhamelRemappin_Result(p).sacc_index              = sacc_index;
        DuhamelRemappin_Result(p).remapping_index         = remapping_index;
        
    end
    
end