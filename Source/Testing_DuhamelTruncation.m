

%
%  Testing_DuhamelTruncation.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_DuhamelTruncation(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-DuhamelTruncation'];
    stimulitype = 'DuhamelTruncation';
    
    % Technical
    dt                              = 0.010; % (s)
    seed                            = 77;
    rng(seed);
    
    % Spatial
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    saccade_threshold               = S_eccentricity/2;
    
    % Temporal
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.300; % (s) w.r.t start of task
    postSaccadefixationPeriod       = 0.300; % (s) time from saccade onset
    %stimuliOffset                  = saccadeOnset;%(s) - turn off stimuli when saccade i
    
    % Utilities - derived
    currentRF                       = 15;%-R_eccentricity:1:R_eccentricity;
    saccades                        = -S_eccentricity:1:S_eccentricity;
    saccadeDelayTime                = roundn((2*S_eccentricity/saccadeSpeed)+0.05,-1); % round to nearest hundred above
    Duration                        = saccadeOnset + saccadeDelayTime + postSaccadefixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    targetOffIntervals{1}           = [];%[stimuliOffset Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(currentRF);
        
        % Rf that will have stim in it, but loose it with saccade
        r = currentRF(i);
        
        % Pick saccade
        s = randi(length(saccades));
        
        % Make sure this saccade is big enough, and has current_ref on retina
        while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= r+saccades(s) && r+saccades(s) <=R_eccentricity))
            s = randi(length(saccades));
        end
        
        % Generate trace
        stimuli{i}.headCenteredTargetLocations  = r;
        stimuli{i}.saccadeTargets               = saccades(s);
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, 0, saccadeOnset, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;
        
        % Add simple information
        stimuli{i}.currentRF = r;
        
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'S_eccentricity', ...
                                    'R_eccentricity', ...
                                    'saccadeOnset', ...
                                    'postSaccadefixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end