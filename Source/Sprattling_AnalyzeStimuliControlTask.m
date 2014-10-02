
%
%  Sprattling_AnalyzeStimuliControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [StimuliControl_Result, Decoded_ReceptiveFieldsLocations, cross_trial_activity] = Sprattling_AnalyzeStimuliControlTask(activity, stimuli)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        %LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-StimuliControl.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-StimuliControl/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;

    % Set parameters
    dt                  = activity.dt;
    R_N                 = activity.R_N;
    R_eccentricity      = stimuli.R_eccentricity;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    stimuliOnset        = stimuli.stimuliOnset;
    stimuliDuration     = stimuli.stimuliDuration;
    
    % Analysis params
    responseWindowStart = 0; % Colby window control, [50ms,250ms] after stim onset.
    responseWindowDuration = 0.200; %(s), colby
    
    latencyWindowSize   = 0.020; % (s), Colby
    latencyWindowLength = timeToTimeStep(latencyWindowSize, dt);

    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Latency & Duration - nonvectorized form is more practical.
    for p=1:numPeriods,
        
        % Find neuron
        neuronIndex = R_eccentricity + stimuli.stimuli{p}.headCenteredTargetLocations + 1;
        
        % Get data for best period of each neuron
        neuronActivity = R_firing_history(neuronIndex, :, p, 1);
        
        % Find latency and duration
        [latencyTimeStep duration] = findNeuronalLatency(neuronActivity, latencyWindowLength);
        
        % Baseline response
        %baseline_response = normalizedIntegration(neuronActivity, dt, 0, stimuliOnset);
        
        % Stim response
        %stim_response   = normalizedIntegration(neuronActivity, dt, stimuliOnset + responseWindowStart, responseWindowDuration);
        
        % Offset response
        %offset_response = normalizedIntegration(neuronActivity, dt, stimuliOnset + stimuliDuration + responseWindowStart, responseWindowDuration);
        
        % Plot
        %figure;plot(neuronActivity);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        
        % Save
        StimuliControl_Result(p).index              = neuronIndex;
        StimuliControl_Result(p).receptiveField     = stimuli.stimuli{p}.headCenteredTargetLocations;
        StimuliControl_Result(p).latency            = stepToTime(latencyTimeStep, dt)-stimuliOnset;
        
        %StimuliControl_Result(p).duration           = duration*dt;
        %StimuliControl_Result(p).baseline_response  = baseline_response;
        %StimuliControl_Result(p).stimulus_response  = stim_response;
        %StimuliControl_Result(p).offset_response    = offset_response;
    end
    
    %% Decode retinal location

    R_Locations = -R_eccentricity:1:R_eccentricity;
    cross_trial_activity = zeros(R_N, numPeriods);
    Decoded_ReceptiveFieldsLocations = zeros(1, R_N);
    
    for i=1:R_N,
        
        % grab all activity for this neuron
        cross_trial_neuronActivity = squeeze(R_firing_history(i, :, :, 1))'; % (R_N, numSavedTimeSteps, numPeriods, numEpochs)
        
        % integrat up activity in each trial
        activity = normalizedIntegration(cross_trial_neuronActivity, dt, stimuliOnset + responseWindowStart, responseWindowDuration); % duration = 200ms , Colby window control, [50ms,250ms] after stim onset.
        
        % decode
        decoded_location = dot(activity, R_Locations)/sum(activity);
        
        % save it
        cross_trial_activity(i, :) = activity;
        Decoded_ReceptiveFieldsLocations(i) = decoded_location;
    end
    
end