
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
        % LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-KusonokiTesting.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-KusonokiTesting/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                          = activity.dt;
    numEpochs                   = activity.numEpochs;
    numPeriods                  = activity.numPeriods;
    
    saccadeOnset                = stimuli.saccadeOnset;
    stimulusOnsetTimes          = stimuli.stimulusOnsetTimes;
    
    % Analysis params
    responseWindowDuration  = 0.300; % (s) from kusonoki paper, it is used in both saccade aligned and stimulus aligned analysis
    stim_responseWindowStart = 0.050 % (s) from kusonoki paper
    
    % Check that this is actually testing stimuli
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Buffers
    stim_buffer = cell(1, length(stimulusOnsetTimes));
    sacc_buffer = cell(1, length(stimulusOnsetTimes));
    
    %% Analysis for each period
    for p=1:numPeriods,
        
        % Rf where stim is located
        stim_location = stimuli.stimuli{k}.headCenteredTargetLocations;
        
        % Get stim onset
        stimOnsetNr = stimuli.stimuli{k}.stimOnsetNr;
        stimuliOnset   = stimuli.stimulusOnsetTimes(stimOnsetNr);
        
        % Get task type, deduce what neuron to record from
        if(stimuli.stimuli{k}.trialType == 1),
            rf = stim_location;
        else
            rf = stim_location - stimuli{k}.saccadeTargets;
        end
        
        % Get neuron index of neuron
        neuronIndex =  R_eccentricity + rf + 1;
        
        % Get data for best period of each neuron
        responseVector  = R_firing_history(neuronIndex, :, p, 1);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(responseVector, dt, stimuliOnset + stim_responseWindowStart, responseWindowDuration);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Save in buffers
        stim_buffer{stimOnsetNr} = [stim_buffer{stimOnsetNr} stimulionset_response];
        sacc_buffer{stimOnsetNr} = [sacc_buffer{stimOnsetNr} saccadeonset_response];

    end
    
    % Turn raw data into struct arrays
    for p=1:numPeriods,
        
        kusonokiSTIMAlignedAnalysis(p).mean = mean(stim_buffer{p});
        kusonokiSTIMAlignedAnalysis(p).std  = std(stim_buffer{p});
        
        kusonokiSACCAlignedAnalysis(p).mean = mean(sacc_buffer{p});
        kusonokiSACCAlignedAnalysis(p).std  = std(sacc_buffer{p});
    end
end