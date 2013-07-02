
%
%  AnalyzeStimuliControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function receptivefield = AnalyzeStimuliControlTask(activityFile, stimuliFile)

    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-StimuliControlTask.mat';
        stimuliFile     = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-StimuliControlTask/stim.mat';
    end
    
    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load(stimuliFile);
    
    % Load data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                  = activity.dt;
    R_N                 = activity.R_N;
    S_N                 = activity.S_N;
    C_N                 = activity.C_N;
    numPeriods          = activity.numPeriods;
    numEpochs           = activity.numEpochs;
    stimuliOnsetDelay   = stimuli.stimuliOnsetDelay;
    onsetTimeStep       = timeToTimeStep(stimuliOnsetDelay, dt);
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Iterate neurons: vectorization beyond this is likely not possible
    
    %% PERIODS, you forgot about that!!!!
    
    % Baseline response
    baseline_activity = R_firing_history(:, 1:onsetTimeStep, :, 1); % [0,saccadeOnset]
    baseline_response = squeeze(trapz(baseline_activity,2));
    
    %should be all identical here for all periods of the same neuron
    baseline_response = baseline_response(:,1); % so we pick the first period for all of them.
    
    % Stimulus response
    stim_activity = R_firing_history(:,timeToTimeStep(stimuliOnsetDelay + 0.050:dt:0.250, dt), :, 1); % [onsetTimeStep+50:250]
    stim_response = squeeze(trapz(stim_activity,2));
    
    % Location
    normalization_response = sum(stim_response,2);
    location = (stim_response*stimuli.headCenteredTargetLocations')./normalization_response';
    
    % Latency
    latency = 2;
    
    % Duration
    
    
    % Make summary figure
    
    % save analysis.mat to directory
    
    receptivefield = 0;
    
end