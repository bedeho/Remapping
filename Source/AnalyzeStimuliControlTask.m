
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
    R_N                 = activity.R_N;
    R_eccentricity      = stimuli.R_eccentricity;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    stimuliOnsetDelay   = stimuli.stimuliOnsetDelay;
    onsetTimeStep       = timeToTimeStep(stimuliOnsetDelay, dt);
    
    % Analysis params
    responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    responseWindowEnd   = responseWindowStart + 0.200;
    
    latencyWindowSize   = 0.020; % (s), Colby
    latencyWindowLength = ceil(latencyWindowSize/dt);
    responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');

    %% Baseline response
    
    %{
    
    % Get time steps in question
    baselineTimeSteps   = 1:(onsetTimeStep-1);
    
    % Extract the given time steps from all neurons in all periods
    baselineActivity    = R_firing_history(:, baselineTimeSteps, :, 1);
    
    % Integrate to find response
    baselineResponse    = squeeze(trapz(baselineActivity,2)); % Integrate
    
    % Should be all identical across neurons and periods, so we pick the first period for all of them.
    baselineResponse    = baselineResponse(:,1);
    
    % Normaliztion step, gives normalized (sp/s) units to response
    baselineResponse    = baselineResponse/(length(baselineTimeSteps)-1);
    
    %% Stimulus response
    
    % Get time steps in question
    activityTimeSteps   = timeToTimeStep(stimuliOnsetDelay + responseWindowStart:dt:responseWindowEnd, dt);
    
    % Extract the given time steps from all neurons in all periods
    stim_activity       = R_firing_history(:, activityTimeSteps, :, 1); % [onsetTimeStep+50:250]
    
    % Integrate to find response
    stim_response       = squeeze(trapz(stim_activity,2));
    
    % Normaliztion step, gives normalized (sp/s) units to response
    stim_response       = stim_response/(length(activityTimeSteps) - 1);
    
    %% Receptive field location
    
    % Level normalization
    stim_response_modenormalized = stim_response - repmat(mode(stim_response, 2),1, numPeriods);
    stim_response_modenormalized(stim_response_modenormalized < 0) = 0;
    
    % Center of mass
    normalization_response = sum(stim_response_modenormalized,2);
    location = (stim_response*stimuli.headCenteredTargetLocations')./normalization_response;
    
    %}
    
    %% Latency & Duration - nonvectorized form is more practical.
    
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
        baselineActivity = neuronActivity(1:(onsetTimeStep-1));
        baselineResponse = squeeze(trapz(baselineActivity,2)); % Integrate
        baselineResponse = baselineResponse/(length(baselineTimeSteps)-1); % Normaliztion step, gives normalized (sp/s) units to response
        
        % Stim response
        activityTimeSteps   = timeToTimeStep(stimuliOnsetDelay + responseWindowStart:dt:responseWindowEnd, dt);  % Get time steps in question
        stim_activity       = neuronActivity(activityTimeSteps); % [onsetTimeStep+50:250]
        stim_response       = squeeze(trapz(stim_activity,2)); % Integrate to find response
        stim_response       = stim_response/(length(activityTimeSteps) - 1); % Normaliztion step, gives normalized (sp/s) units to response
        
        % Plot
        figure;plot(neuronActivity);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
        
        % Save
        StimuliControl_Neurons(c).index            = neuronIndex;
        StimuliControl_Neurons(c).latencyTimeStep  = stepToTime(latencyTimeStep, dt);
        StimuliControl_Neurons(c).duration         = duration*dt;
        StimuliControl_Neurons(c).baselineResponse = baselineResponse;
        StimuliControl_Neurons(c).stimulusresponse = stim_response;
        
        StimuliControl_indexes(c) = neuronIndex;

        c = c + 1;
    end
    
    x=1;
end