
%
%  AnalyzeDuhamelTruncation.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelTruncation_Neurons, DuhamelTruncation_indexes] = AnalyzeDuhamelTruncation(activity, stimuli)

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
    
    % Analysis params
    %latencyWindowSize   = 0.020; % (s), colby papers
    %latencyWindowLength = ceil(latencyWindowSize/dt);
    RF_inclusion_th     = 5; % (deg) neurons this far away from any given trial are analysed together
    responseWindowDuration = 0.200;
    %responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    
    c = 1;
    DuhamelTruncation_indexes = [];
    
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation_HeadCentered = stimuli.headCenteredTargetLocations(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,currentRFLocation_HeadCentered - RF_inclusion_th):1:min(R_eccentricity,currentRFLocation_HeadCentered + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data neuron
            neuronActivity  = R_firing_history(neuronIndex, :, p, 1);
            
            % Offset response
            saccadeonset_response = normalizedIntegration(neuronActivity, dt, saccadeOnset, responseWindowDuration);
            
            %% DEBUG - looked good
            %figure;plot(neuronActivity);hold on; plot(timeToTimeStep([saccadeOnset saccadeOnset], dt), [0 1], 'r');
            
            % Save
            DuhamelTruncation_Neurons(c).index                   = neuronIndex;
            DuhamelTruncation_Neurons(c).saccadeonset_response   = saccadeonset_response;
            
            DuhamelTruncation_indexes(c) = neuronIndex;
            
            c = c + 1;
            
        end
    end
    
end