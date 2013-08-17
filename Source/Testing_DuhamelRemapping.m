
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
        stimuliOnset                    = saccadeOnset - stimuliDuration; % (s) w.r.t start of task
        postSaccadefixationPeriod       = 0.300; % (s) time from saccade COMPLETION ESTIMATE "saccadeDelayTime" below.
    end
    
    %% Utilities - derived
    futureRepetiveField             = -R_eccentricity:1:R_eccentricity; % Remapping TARGET, i.e. postsaccadic (-R_eccentricity+R_edge_effect_buffer):1:(R_eccentricity+R_edge_effect_buffer)
    saccades                        = -S_eccentricity:1:S_eccentricity; % Pick among these saccades
    saccadeDelayTime                = roundn((2*S_eccentricity/saccadeSpeed) + 0.05,-1); % round to nearest hundred above
    Duration                        = saccadeOnset + saccadeDelayTime + postSaccadefixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    targetOffIntervals              = {[0 stimuliOnset]; [(stimuliOnset+stimuliDuration) Duration]}
    
    %% Generate stimuli
    for i = 1:length(futureRepetiveField);
        
        % Location of target
        r = futureRepetiveField(i);
        
        % Pick saccade
        s = randi(length(saccades));
        
        % Make sure it is big enough and keeps target on retina
        while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= r-saccades(s) && r-saccades(s) <=R_eccentricity)) % Make sure this saccade is big enough, and has current_ref on retina
            s = randi(length(saccades));
        end
        
        stimuli{i}.headCenteredTargetLocations  = r-saccades(s);
        stimuli{i}.saccadeTargets               = saccades(s);
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.numSaccades                  = length(saccadeOnset);
        
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, 0, saccadeOnset, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;

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
                                    'fixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end