
%
%  AnalyzeStimuliControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [StimuliControl_Neurons, StimuliControl_indexes] = AnalyzeStimuliControlTask(activity, stimuli)

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
    R_eccentricity      = stimuli.R_eccentricity;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    stimuliOnsetDelay   = stimuli.stimuliOnsetDelay;
    stimuliDuration     = stimuli.stimuliDuration;
    
    % Analysis params
    responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    responseWindowDuration = 0.200; %(s), colby
    %responseWindowEnd   = responseWindowStart + responseWindowDuration;
    
    latencyWindowSize   = 0.020; % (s), Colby
    latencyWindowLength = ceil(latencyWindowSize/dt);
    responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Latency & Duration - nonvectorized form is more practical.
    
    c = 1;
    StimuliControl_indexes = [];
    
    for p=1:numPeriods,
        
        % Find neuron
        neuronIndex = R_eccentricity + stimuli.headCenteredTargetLocations(p) + 1;
        
        % Get data for best period of each neuron
        neuronActivity = R_firing_history(neuronIndex, :, p, 1);
        
        % Find latency and duration
        [latencyTimeStep duration] = findNeuronalLatency(responseThreshold, neuronActivity, latencyWindowLength);
        
        % Baseline response
        baseline_response = normalizedIntegration(neuronActivity, dt, 0, stimuliOnsetDelay);
        
        % Stim response
        stim_response   = normalizedIntegration(neuronActivity, dt, stimuliOnsetDelay + responseWindowStart, responseWindowDuration);
        
        % Offset response
        offset_response = normalizedIntegration(neuronActivity, dt, stimuliOnsetDelay + stimuliDuration + responseWindowStart, responseWindowDuration);
        
        % Plot
        %figure;plot(neuronActivity);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        
        % Save
        StimuliControl_Neurons(c).index            = neuronIndex;
        StimuliControl_Neurons(c).latency          = stepToTime(latencyTimeStep, dt)-stimuliOnsetDelay;
        StimuliControl_Neurons(c).duration         = duration*dt;
        StimuliControl_Neurons(c).baselineResponse = baseline_response;
        StimuliControl_Neurons(c).stimulusresponse = stim_response;
        StimuliControl_Neurons(c).offset_response  = offset_response;
        
        StimuliControl_indexes(c) = neuronIndex;

        c = c + 1;
    end
end