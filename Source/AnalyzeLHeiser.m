
%
%  AnalyzeLHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 13/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function [LHeiser_Result] = AnalyzeLHeiser(LHeiser_DuhamelRemappingTrace, stimuli)
    
    % Get data remapping activity data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    S_eccentricity       = stimuli.S_eccentricity;
    saccadeOnset         = stimuli.saccadeOnset;
    stimuliOnset         = stimuli.stimuliOnset;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    
    stim_stimuliOnset = stim_stimuli.stimuliOnset;
    sacc_saccadeOnset = sacc_stimuli.saccadeOnset;
    
    % Analysis params
    responseWindowDuration = 0.200;
    
    stim_responseWindowStart = 0.050; % Colby window control, [50ms,250ms] after stim onset.
    
    % For plots
    stimuliOnsetTimeStep = timeToTimeStep(stimuliOnset, dt);
    saccadeOnsetTimeStep = timeToTimeStep(saccadeOnset, dt);
    stim_stimuliOnsetTimeStep = timeToTimeStep(stim_stimuliOnset, dt);
    sacc_saccadeOnsetTimeStep = timeToTimeStep(sacc_saccadeOnset, dt);
    
    responseWindowSize = timeToTimeStep(responseWindowDuration, dt);
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(size(stim_control_activity,3) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');
    assert(size(sacc_control_activity,3) == 2*S_eccentricity+1, 'Saccade control must be done for saccade map');
    
    %% Latency & Duration
    for p=1:numPeriods,
        
        % Index of R neuron which will recieve remapping activity: currentRF
        remappedInto_neuronIndex = R_eccentricity + stimuli.stimuli{p}.currentRF + 1;
        
        % Get data for neuron
        remap_responseVector  = R_firing_history(remappedInto_neuronIndex, :, p, 1);

        % Find latency and duration
        latencyTimeStep = findNeuronalLatency_NEW(remap_responseVector, dt);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(remap_responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(remap_responseVector, dt, stimuliOnset, responseWindowDuration);
        
        %% Save
        LHeiser_Result(p).index                   = remappedInto_neuronIndex;
        LHeiser_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        LHeiser_Result(p).futureRF                = stimuli.stimuli{p}.futureRF;
        LHeiser_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        LHeiser_Result(p).remappingLatency        = (latencyTimeStep - timeToTimeStep(saccadeOnset,dt))*dt;
        LHeiser_Result(p).stimLatency             = (stim_latencyTimeStep - timeToTimeStep(stim_stimuliOnset, dt))*dt;
        
        LHeiser_Result(p).saccadeonset_response   = saccadeonset_response;
        LHeiser_Result(p).stimulionset_response   = stimulionset_response;
        
        LHeiser_Result(p).stim_index              = stim_index;
        LHeiser_Result(p).sacc_index              = sacc_index;
        LHeiser_Result(p).remapping_index         = remapping_index;
        
    end
    
end