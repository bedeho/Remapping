
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
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    saccadeOnset         = stimuli.saccadeOnset;
    stimuliOnset         = stimuli.stimuliOnset;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    
    % Analysis params
    latencyWindowSize   = 0.020; % (s), colby papers
    latencyWindowLength = timeToTimeStep(latencyWindowSize, dt);
    responseWindowDuration = 0.200;
    responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation_HeadCentered = stimuli.stimuli{p}.headCenteredTargetLocations;
        futureRFLocation = currentRFLocation_HeadCentered - stimuli.stimuli{p}.saccadeTargets;
        
        % Find neuron
        neuronIndex = R_eccentricity + futureRFLocation + 1;

        % Get data for best period of each neuron
        responseVector  = R_firing_history(neuronIndex, :, p, 1);

        % Find latency and duration
        [latencyTimeStep, duration] = findNeuronalLatency(responseThreshold, responseVector, latencyWindowLength);
        
        %{
        figure;
        subplot(1,2,1);
        plot(responseVector);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        subplot(1,2,2);
        plot(responseVector);hold on; imagesc(R_firing_history(:, :, p, 1));
        %}
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(responseVector, dt, stimuliOnset, responseWindowDuration);
        
        
        stim_control_activity, sacc_control_activity
        
        
        % Save
        DuhamelRemappin_Result(p).index                   = neuronIndex;
        DuhamelRemappin_Result(p).futureRF                = futureRFLocation;
        DuhamelRemappin_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        
        DuhamelRemappin_Result(p).latency                 = stepToTime(latencyTimeStep, dt)-saccadeOnset;
        DuhamelRemappin_Result(p).Duration                = duration*dt;
        DuhamelRemappin_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelRemappin_Result(p).stimulionset_response   = stimulionset_response;
    end
    
end