
%
%  Sprattling_AnalyzeDuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelRemappin_Result, Decoded_ReceptiveFieldsLocations] = Sprattling_AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli, Decoded_ReceptiveFieldsLocations)

    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
       
        %LoadActivity
        exp = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/sprattling_visual_learning/';
        
        activity = load([exp 'baseline/TrainedNetwork/activity-basic-DuhamelRemappingTrace.mat']);
        stimuli  = load([exp 'STIM-basic-DuhamelRemappingTrace/stim.mat']);
        
        x = load([exp 'baseline/TrainedNetwork/activity-basic-StimuliControl.mat']);
        stim_control_activity = x.R_firing_history;
        stim_stimuli = load([exp 'STIM-basic-StimuliControl/stim.mat']);
        
        x = load([exp 'baseline/TrainedNetwork/activity-basic-SaccadeControl.mat'])
        sacc_control_activity = x.R_firing_history;
        sacc_stimuli = load([exp 'STIM-basic-SaccadeControl/stim.mat']);
        
    end
    
    % Get data remapping activity data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_N                  = activity.R_N;
    R_eccentricity       = stimuli.R_eccentricity;
    S_eccentricity       = stimuli.S_eccentricity;
    saccadeOnset         = stimuli.saccadeOnset;
    stimuliOnset         = stimuli.stimuliOnset;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    
    stim_stimuliOnset = stim_stimuli.stimuliOnset;
    sacc_saccadeOnset = sacc_stimuli.saccadeOnset;
    
    % Analysis params
    responseWindowDuration = 0.300; % CLASSIC:0.200, but LHeiser2005 uses 0.300
    
    % Check that we have done the right controls and that this is truly a
    % testing data set (epoch==1)
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    assert(size(stim_control_activity,3) == 2*R_eccentricity+1, 'Stimulus control must be done for entire retina');
    assert(size(sacc_control_activity,3) == 2*S_eccentricity+1, 'Saccade control must be done for saccade map');
    
    % Latency & Duration
    for p=1:numPeriods,
        
        % Index of R neuron which will recieve remapping activity: currentRF
        [remappedInto_neuronIndex, loc] = get_neuron_from_RF(stimuli.stimuli{p}.currentRF); % R_eccentricity + stimuli.stimuli{p}.currentRF + 1;
        
        % Get data for neuron
        remap_responseVector = R_firing_history(remappedInto_neuronIndex, :, p, 1);

        % Find latency and duration
        %[latencyTimeStep, duration] = findNeuronalLatency(remap_responseVector, latencyWindowLength);
        latencyTimeStep = findNeuronalLatency_NEW(remap_responseVector, dt);
        
        % Saccade aligned response window
        saccadeonset_response = normalizedIntegration(remap_responseVector, dt, saccadeOnset, responseWindowDuration);
        
        % Stimuli aligned response window
        stimulionset_response = normalizedIntegration(remap_responseVector, dt, stimuliOnset, responseWindowDuration);
        
        %% STIM - CONTROL
        
        % Stimuli Control Trial number where stimuli is presented at stimuli.stimuli{p}.futureRF
        futureRF_neuronIndex = R_eccentricity + stimuli.stimuli{p}.futureRF + 1;
        
        % Stim response of remapped neuron when stim is in future RF
        stim_responseVector = stim_control_activity(remappedInto_neuronIndex, :, futureRF_neuronIndex);
        stim_responseWindowStart = stim_stimuliOnset + 0.200; % CLASIC: +0.050, but we use LHeiser2005: 0.200
        stim_control_response = normalizedIntegration(stim_responseVector, dt, stim_responseWindowStart, responseWindowDuration);
        stim_index = saccadeonset_response - stim_control_response;
        
        % Find latency and duration
        stim2_responseVector = stim_control_activity(futureRF_neuronIndex, :, futureRF_neuronIndex);
        %[stim_latencyTimeStep stim_duration] = findNeuronalLatency(stim2_responseVector, latencyWindowLength);
        stim_latencyTimeStep= findNeuronalLatency_NEW(stim2_responseVector, dt);
        
        %% SACC = CONTROL
        
        % Index of R neuron which has RF representing the saccade
        % performed
        sacc_neuronIndex = S_eccentricity + stimuli.stimuli{p}.saccadeTargets + 1;
        
        % Sacc response of remapped neuron when given saccade is performed
        sacc_responseVector = sacc_control_activity(remappedInto_neuronIndex, :, sacc_neuronIndex);
        saccade_control_response = normalizedIntegration(sacc_responseVector, dt, sacc_saccadeOnset, responseWindowDuration);
        sacc_index = saccadeonset_response - saccade_control_response;
        
        % Remapping Indexes
        if(stim_index > 0 && sacc_index > 0),
            remapping_index = sqrt(stim_index^2 + sacc_index^2);
            %remapping_index = (stim_index + sacc_index)/2;
        else
            remapping_index = 0;
        end
        
        %% Save
        DuhamelRemappin_Result(p).index                   = remappedInto_neuronIndex;
        DuhamelRemappin_Result(p).currentRF               = stimuli.stimuli{p}.currentRF;
        DuhamelRemappin_Result(p).futureRF                = stimuli.stimuli{p}.futureRF;
        DuhamelRemappin_Result(p).saccade                 = stimuli.stimuli{p}.saccadeTargets;
        DuhamelRemappin_Result(p).remappingLatency        = (latencyTimeStep - timeToTimeStep(saccadeOnset,dt))*dt;
        DuhamelRemappin_Result(p).stimLatency             = (stim_latencyTimeStep - timeToTimeStep(stim_stimuliOnset, dt))*dt;
        
        %DuhamelRemappin_Result(p).Duration                = duration*dt;
        DuhamelRemappin_Result(p).saccadeonset_response   = saccadeonset_response;
        DuhamelRemappin_Result(p).stimulionset_response   = stimulionset_response;
        
        DuhamelRemappin_Result(p).stim_index              = stim_index;
        DuhamelRemappin_Result(p).sacc_index              = sacc_index;
        DuhamelRemappin_Result(p).remapping_index         = remapping_index;
        
    end
    
    function [idx, loc] = get_neuron_from_RF(location)
        
        % find closest
        [abs_diff, idx] = min(abs(Decoded_ReceptiveFieldsLocations-location));
        
        % get location
        loc = Decoded_ReceptiveFieldsLocations(idx);
        
        assert(abs(location-loc) <= 5, 'Not possible to find neuron close enough');
        %disp('Not possible to find neuron close enough');
        
    end
    
end