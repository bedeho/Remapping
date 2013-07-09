
%
%  AnalyzeDuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [NeuronsAnalyzed, LatencyTimeSteps, Durations] = AnalyzeDuhamelRemapping(activity, stimuli)

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
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    
    % Analysis params
    latencyWindowSize   = 0.020; % (s), colby papers
    latencyWindowLength = ceil(latencyWindowSize/dt);
    RF_inclusion_th     = 5; % (deg) neurons this far away from any given trial are analysed together
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    
    c = 0;
    NeuronsAnalyzed = [];
    LatencyTimeSteps = [];
    Durations = [];
    
    % Allocate some result buffers
    %foundOnset       = zeros(1, R_N);
    %foundOffset      = zeros(1, R_N);
    %latencyTimeStep  = nan*zeros(1, R_N);
    %durationTimeStep = nan*zeros(1, R_N);
    
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation = stimuli.headCenteredTargetLocations(p);
        futureRFLocation = currentRFLocation + stimuli.saccadeTargets(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,futureRFLocation - RF_inclusion_th):1:min(R_eccentricity,futureRFLocation + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data for best period of each neuron
            responseVector  = R_firing_history(neuronIndex, :, p, 1);

            % Find start of vector
            [latencyTimeStep, duration] = findNeuronalLatency(responseThreshold, responseVector, latencyWindowLength);
            
            % Save
            NeuronsAnalyzed(c)  = neuronIndex;
            LatencyTimeSteps(c) = timeStepToTime(latencyTimeStep,dt);
            Duration(c)         = duration*dt;
            
            c = c + 1;
            
        end
    end
    
end