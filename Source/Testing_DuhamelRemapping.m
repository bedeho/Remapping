
%
%  Testing_DuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 06/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_DuhamelRemapping(Name, stimulitype, saccadeOnset, stimuliDuration, stimuliOnset, postSaccadefixationPeriod)

    % Import global variables
    declareGlobalVars();
    global base;
    
    %% Parameters
    
    % Technical
    dt                              = 0.010; % (s)
    seed                            = 77;
    rng(seed);
    
    % Spatial
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    saccade_threshold               = S_eccentricity/2;
    
    % Dynamical
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    
    % Temporal
    if nargin<6,
        stimulitype                     = 'DuhamelRemapping';
        saccadeOnset                    = 0.300; % (s) w.r.t start of task
        stimuliDuration                 = 0.100; % (s)
        stimuliOnset                    = saccadeOnset - 2*stimuliDuration; % (s) w.r.t start of task
        postSaccadefixationPeriod       = 0.300; % (s) time from saccade COMPLETION ESTIMATE "saccadeDelayTime" below.
    end
    
    %% Utilities - derived
    currentRF  = 10;%-R_eccentricity:1:R_eccentricity; % Remapping TARGET, i.e. postsaccadic (-R_eccentricity+R_edge_effect_buffer):1:(R_eccentricity+R_edge_effect_buffer)
    saccades                        = -S_eccentricity:1:S_eccentricity; % Pick among these saccades
    saccadeDelayTime                = roundn((2*S_eccentricity/saccadeSpeed) + 0.05,-1); % round to nearest hundred above
    Duration                        = saccadeOnset + saccadeDelayTime + postSaccadefixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    targetOffIntervals{1}           = [0 stimuliOnset; (stimuliOnset+stimuliDuration) Duration];
    
    %% Generate stimuli
    for i = 1:length(currentRF);
        
        % Rf that will get remapping activity INTO it
        r = currentRF(i);
        
        % Pick saccade
        s = randi(length(saccades));
        
        % Make sure it is big enough and keeps target on retina
        while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= r+saccades(s) && r+saccades(s) <=R_eccentricity)) % Make sure this saccade is big enough, and has current_ref on retina
            s = randi(length(saccades));
        end
        
        % Screen location of stimulus
        stim_screen_location = r+saccades(s);
        
        % Generate trace
        stimuli{i}.headCenteredTargetLocations  = stim_screen_location;
        stimuli{i}.saccadeTargets               = saccades(s);
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, stimuli{i}.targetOffIntervals, 0, stimuli{i}.saccadeTimes, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;
        
        % Add simple information
        stimuli{i}.currentRF             = r;
        stimuli{i}.stim_screen_location  = stim_screen_location;
        stimuli{i}.numSaccades           = length(saccadeOnset);
    end
    
    % Save params
    filename = [Name '-' stimulitype];
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'S_eccentricity', ...
                                    'R_eccentricity', ...
                                    'saccadeOnset', ...
                                    'stimuliDuration', ...
                                    'stimuliOnset', ...
                                    'postSaccadefixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end