
%
%  AnalyzeStimuliControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [baselineResponse, stim_response, location, foundOnset, foundOffset, latencyTimeStep, durationTimeStep, neuronResponse] = AnalyzeStimuliControlTask(activity, stimuli)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        activity = LoadActivity('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-StimuliControlTask.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-StimuliControlTask/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                  = activity.dt;
    R_N                 = activity.R_N;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    stimuliOnsetDelay   = stimuli.stimuliOnsetDelay;
    onsetTimeStep       = timeToTimeStep(stimuliOnsetDelay, dt);
    
    % Analysis params
    responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    responseWindowEnd   = responseWindowStart + 0.200;
    latencyWindowSize   = 0.020; % (s),
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');

    %% Baseline response
    
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
    
    %% Latency & Duration
    
    % vectorization became problematic when it came to extracting the right
    % periods and neurons simultanously, so I skipped it.
    
    % Allocate some result buffers
    foundOnset       = zeros(1, R_N);
    foundOffset      = zeros(1, R_N);
    latencyTimeStep  = nan*zeros(1, R_N);
    durationTimeStep = nan*zeros(1, R_N);
    neuronResponse   = zeros(R_N, size(R_firing_history,2));
    
    for n=1:R_N,
        
        dist = abs(location(n) - stimuli.headCenteredTargetLocations);
        [C,I] = min(dist);
        
        % Get data for best period of each neuron
        bestPeriodActivity  = R_firing_history(n, :, I, 1);
        
        % Save to figure
        neuronResponse(n, :) = bestPeriodActivity;
        
        % Compute response in different latency shifts
        thresholdResponse   = stim_response(n, I)*0.5; % stimulus period when stimuli is closest
        
        num = length(bestPeriodActivity);
        windowResponse = zeros(1, num);
        latencyWindowOffset = 0:ceil(latencyWindowSize/dt);

        for t=onsetTimeStep:(length(bestPeriodActivity)-length(latencyWindowOffset)+1),
            
            % Get time steps in question
            latencyWindow = t + latencyWindowOffset;
            
            % Integrate to find response
            windowResponse(t) = trapz(bestPeriodActivity(latencyWindow));
            
            % Normalize
            windowResponse(t) = windowResponse(t)/(length(latencyWindow)-1);
            
            if ~foundOnset(n),
                
                if windowResponse(t) >= thresholdResponse;
                    latencyTimeStep(n) = t;
                    foundOnset(n) = true;
                end
            elseif ~foundOffset(n),
                
                if windowResponse(t) <= thresholdResponse;
                    durationTimeStep(n) = t;
                    foundOffset(n) = true;
                end
            end
        end        
    end
    
end