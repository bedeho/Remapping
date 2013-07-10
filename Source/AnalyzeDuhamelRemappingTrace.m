
%
%  AnalyzeDuhamelRemappingTrace.m
%  Remapping
%
%  Created by Bedeho Mender on 10/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelRemapping_analyzedNeurons, DuhamelRemapping_indexes] = AnalyzeDuhamelRemappingTrace(activity, stimuli)

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
    responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    
    c = 1;
    DuhamelRemapping_neuronIndexes = [];
    
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation_HeadCentered = stimuli.headCenteredTargetLocations(p);
        futureRFLocation = currentRFLocation_HeadCentered - stimuli.saccadeTargets(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,futureRFLocation - RF_inclusion_th):1:min(R_eccentricity,futureRFLocation + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data for best period of each neuron
            responseVector  = R_firing_history(neuronIndex, :, p, 1);
            
            % Find latency and duration
            [latencyTimeStep, duration] = findNeuronalLatency(responseThreshold, responseVector, latencyWindowLength);
            
            %figure;plot(responseVector);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
            
            % Save
            DuhamelRemapping_Neurons(c).index            = neuronIndex;
            DuhamelRemapping_Neurons(c).latencyTimeStep  = stepToTime(latencyTimeStep, dt);
            DuhamelRemapping_Neurons(c).Duration         = duration*dt;
            
            DuhamelRemapping_indexes(c) = neuronIndex;
            
            c = c + 1;
            
        end
    end
    
end