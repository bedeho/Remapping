
%
%  AnalyzeStimuliControlTask.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function receptivefield = AnalyzeStimuliControlTask(activityFile, stimuliFile)

    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/-R_w_INHB=0.10989/PrewiredNetwork/activitykusonoki.mat';
        stimuliFile     = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-KusonokiTesting/stim.mat';
    end
    
    %{
    
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
    stimuliOnsetDelay   = stimuli.saccadeOnset;
    onsetTimeStep       = timeToTimeStep(stimuliOnsetDelay, dt);
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% PERIODS, you forgot about that!!!!
    
    % Baseline response
    baseline_activity = R_firing_history(:, 0:onsetTimeStep, :, 1); % [0,saccadeOnset]
    baseline_response = trapz(baseline_activity')';
    
    % Stimulus response
    stim_activity = R_firing_history(:, onsetTimeStep + dt*(50:250), :, 1); % [onsetTimeStep+50:250]
    stim_response = trapz(baseline_activity')';
    
    % Latency
    
    % Duration
    
    % Location
    
    % Make summary figure
    
    % save analysis.mat to directory
    
    %}
    
    receptivefield = 0;
    
end