
%
%  AnalyzeSaccadeControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 01/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [SaccadeControl_Result] = AnalyzeSaccadeControlTask(activity, stimuli)

    % Check if this is manual run
    if nargin == 0,
        
        disp('Loading input files...');
        % LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-SaccadeControlTask.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-SaccadeControlTask/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                  = activity.dt;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    saccadeOnset        = stimuli.saccadeOnset;
    S_eccentricity      = stimuli.S_eccentricity;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
        
    % Analysis params
    responseWindowDuration  = 0.200; % (s), Colby window aligned at saccadeOFFSET
    
    % Turn into struct array
    for p=1:numPeriods,
        
        saccade = stimuli.stimuli{p}.saccadeTargets;
        index   = S_eccentricity + saccade + 1;
        
        SaccadeControl_Result(p).index            = index;
        SaccadeControl_Result(p).receptiveField   = saccade;
        SaccadeControl_Result(p).saccadeonset_response = normalizedIntegration(R_firing_history(index,:,p), dt, saccadeOnset, responseWindowDuration); % [onsetTimeStep+50:250]
    end
end