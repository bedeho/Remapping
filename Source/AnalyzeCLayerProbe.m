
%
%  AnalyzeCLayerProbe.m
%  Remapping
%
%  Created by Bedeho Mender on 13/08/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity, stimuli)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        %LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-CLayerProbe.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-CLayerProbe/stim.mat');
    end
    
    % Get data
    C_firing_history_flat= activity.C_firing_history_flat;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    C_N                  = activity.C_N;
    Duration             = stimuli.Duration;
    
    headCenteredTargetLocations = stimuli.headCenteredTargetLocations;
    saccadeTargets              = stimuli.saccadeTargets;
    
    numVisualLocations   = length(headCenteredTargetLocations);
    numSaccades          = length(saccadeTargets);
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Record responses
    C_responses = zeros(numSaccades, numVisualLocations, C_N);
    
    p = 1;
    for i = 1:numVisualLocations,
            
        for j = 1:numSaccades,
            
            % Layer activity
            %neuronActivity  = C_firing_history(:, :, p, 1);
            
            % Compute and save response
            C_responses(j, i, :) =  C_firing_history_flat(:, p); %normalizedIntegration(neuronActivity, dt, 0, Duration);

            % Next period
            p = p + 1;
        end
    end
    
    %% Compute preferences
    
    CLabeProbe_Neurons_S = zeros(1, C_N);
    CLabeProbe_Neurons_V = zeros(1, C_N);
    
    [targetMesh saccadeMesh] = meshgrid(headCenteredTargetLocations, saccadeTargets);
    
    % Could have vectorized outer loop, but who cares
    for c=1:C_N,
        
        % Get neuron response
        response = C_responses(:, :, c);
        
        % Find center of mass
        norm = sum(sum(response));
        
        CLabeProbe_Neurons_S(c) = sum(sum(response .* saccadeMesh))/norm;
        CLabeProbe_Neurons_V(c) = sum(sum(response .* targetMesh))/norm;
        
    end
    
end