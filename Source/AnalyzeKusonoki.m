
%
%  AnalyzeKusonoki.m
%  Remapping
%
%  Created by Bedeho Mender on 01/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activity, stimuli)

    % Check if this is manual run
    if nargin == 0,
        
        disp('Loading input files...');
        activity = LoadActivity('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-KusonokiTesting.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-KusonokiTesting/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                          = activity.dt;
    R_N                         = activity.R_N;
    numEpochs                   = activity.numEpochs;
    numPeriods                  = activity.numPeriods;
    
    saccadeOnset                = stimuli.saccadeOnset;
    stimulusOnsetTimes          = stimuli.stimulusOnsetTimes;
    saccadeTargets              = stimuli.saccadeTargets;
    headCenteredTargetLocations = stimuli.headCenteredTargetLocations;
    
    numStimuliOnsetTimes        = length(stimulusOnsetTimes);
    numSaccadeTargets           = length(saccadeTargets);
    numTargetsLocations         = length(headCenteredTargetLocations);
                  
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Analysis params
    responseWindowSize  = 0.300; % (s) from kusonoki paper, it is used in both saccade aligned and stimulus aligned analysis
    
    %% Analysis for each period
    for p=1:numPeriods,
        
        % Params of this period
        targetNr    = stimuli.stimuli{p}.targetNr;
        saccadeNr   = stimuli.stimuli{p}.saccadeNr;
        stimOnsetNr = stimuli.stimuli{p}.stimOnsetNr;
        
        stimOnsetTime        = stimulusOnsetTimes(stimOnsetNr);
        saccOnsetTimeSteps   = timeToTimeStep(saccadeOnset + 0:dt:responseWindowSize, dt);
        stimOnsetTimeSteps   = timeToTimeStep(stimOnsetTime + 0:dt:responseWindowSize, dt);
        
        % Extract the given time steps from all neurons in all periods
        saccade_activity    = R_firing_history(:, saccOnsetTimeSteps, :, 1);
        stimulus_activity   = R_firing_history(:, stimOnsetTimeSteps, :, 1);
        
        % Integrate to find response
        saccade_response    = squeeze(trapz(saccade_activity,2));
        stimulus_response   = squeeze(trapz(stimulus_activity,2));

        % Normaliztion step, gives normalized (sp/s) units to response
        saccade_response    = saccade_response/(length(saccOnsetTimeSteps) - 1);
        stimulus_response   = stimulus_response/(length(stimOnsetTimeSteps) - 1);            
        
        % Save results
        kusonokiSACCAlignedAnalysis(p).targetNr     = targetNr;
        kusonokiSACCAlignedAnalysis(p).saccadeNr    = saccadeNr;
        kusonokiSACCAlignedAnalysis(p).stimOnsetNr  = stimOnsetNr;
        kusonokiSACCAlignedAnalysis(p).response     = saccade_response;
        
        kusonokiSTIMAlignedAnalysis(p).targetNr     = targetNr;
        kusonokiSTIMAlignedAnalysis(p).saccadeNr    = saccadeNr;
        kusonokiSTIMAlignedAnalysis(p).stimOnsetNr  = stimOnsetNr;
        kusonokiSTIMAlignedAnalysis(p).response     = stimulus_response;
        
    end
end