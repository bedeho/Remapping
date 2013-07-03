
%
%  AnalyzeSaccadeControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 01/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function saccade_response = AnalyzeSaccadeControlTask(activityFile, stimuliFile)

    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-SaccadeControlTask.mat';
        stimuliFile     = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-SaccadeControlTask/stim.mat';
    end
    
    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load(stimuliFile);
    
    % Load data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                  = activity.dt;
    numEpochs           = activity.numEpochs;
    saccadeOnsetDelay   = stimuli.saccadeOnsetDelay;

    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
        
    % Analysis params
    responseWindowSize  = 0.200; % (s), Colby window aligned at saccadeOFFSET
    
    %% Saccade onset
    
    % Get time steps in question
    activityTimeSteps   = timeToTimeStep(saccadeOnsetDelay + 0:dt:responseWindowSize, dt);
    
    % Extract the given time steps from all neurons in all periods
    saccade_activity    = R_firing_history(:, activityTimeSteps, :, 1); % [onsetTimeStep+50:250]
    
    % Integrate to find response
    saccade_response    = squeeze(trapz(saccade_activity,2));
    
    % Normaliztion step, gives normalized (sp/s) units to response
    saccade_response    = saccade_response/(length(activityTimeSteps) - 1);
    
end