
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
    R_eccentricity              = stimuli.R_eccentricity;
    numEpochs                   = activity.numEpochs;
    numPeriods                  = activity.numPeriods;
    
    saccadeOnset                = stimuli.saccadeOnset;
    stimulusOnsetTimes          = stimuli.stimulusOnsetTimes;
    
    % Analysis params
    responseWindowDuration  = 0.300; % (s) from kusonoki paper, it is used in both saccade aligned and stimulus aligned analysis
    stim_responseWindowStart = 0.050; % (s) from kusonoki paper
    
    % Check that this is actually testing stimuli
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    % Buffers
    stim_buffer = cell(2, length(stimulusOnsetTimes));
    sacc_buffer = cell(2, length(stimulusOnsetTimes));
    
    %% Analysis for each period
    for p=1:numPeriods,
        
        % Rf where stim is located
        stim_location = stimuli.stimuli{p}.headCenteredTargetLocations;
        
        % Get stim onset
        stimOnsetNr  = stimuli.stimuli{p}.stimOnsetNr;
        stimuliOnset = stimuli.stimulusOnsetTimes(stimOnsetNr);
        
        % Get task type, deduce what neuron to record from
        if(stimuli.stimuli{p}.trialType == 1),
            rf = stim_location;
        else
            rf = stim_location - stimuli.stimuli{p}.saccadeTargets;
        end
        
        % Get neuron index of neuron
        neuronIndex =  R_eccentricity + rf + 1;
        
        % Get data for best period of each neuron
        responseVector  = R_firing_history(neuronIndex, :, p, 1);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(responseVector, dt, stimuliOnset + stim_responseWindowStart, responseWindowDuration);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Get task type, and save appripriately
        if(stimuli.stimuli{p}.trialType == 1),
            stim_buffer{1,stimOnsetNr} = [stim_buffer{1,stimOnsetNr} stimulionset_response];
            sacc_buffer{1,stimOnsetNr} = [sacc_buffer{1,stimOnsetNr} saccadeonset_response];
        else
            stim_buffer{2,stimOnsetNr} = [stim_buffer{2,stimOnsetNr} stimulionset_response];
            sacc_buffer{2,stimOnsetNr} = [sacc_buffer{2,stimOnsetNr} saccadeonset_response];
        end
    end
    
    % Turn raw data into struct arrays
    for i=1:length(stimulusOnsetTimes),
        
        % stim aligned analysis
        kusonokiSTIMAlignedAnalysis(i).current_mean = mean(stim_buffer{1,i});
        kusonokiSTIMAlignedAnalysis(i).current_std  = std(stim_buffer{1,i});
        
        kusonokiSTIMAlignedAnalysis(i).future_mean = mean(stim_buffer{2,i});
        kusonokiSTIMAlignedAnalysis(i).future_std  = std(stim_buffer{2,i});
        
        % sacca aligned analysis
        kusonokiSACCAlignedAnalysis(i).current_mean = mean(sacc_buffer{1,i});
        kusonokiSACCAlignedAnalysis(i).current_std  = std(sacc_buffer{1,i});
        
        kusonokiSACCAlignedAnalysis(i).future_mean = mean(sacc_buffer{2,i});
        kusonokiSACCAlignedAnalysis(i).future_std  = std(sacc_buffer{2,i});
    end
end