

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
    
    % Params
    dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    
    saccade_threshold               = S_eccentricity/2;
    
    % Dynamical quantities
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.300; % (s) w.r.t start of task
    fixationPeriod                  = 0.300; % (s) time from saccade onset
    stimuliOffset                   = saccadeOnset/2;%0.050; % (s)
    
    % Generate stimuli
    rng(seed);
    saccadeDelayTime                = roundn((2*S_eccentricity/saccadeSpeed)+0.05,-1) % round to nearest hundred above
    Duration                        = saccadeOnset + saccadeDelayTime + fixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    headCenteredTargetLocations     = 0;%-R_eccentricity:1:R_eccentricity;
    saccadeTargets                  = zeros(1, length(headCenteredTargetLocations));

    
    saccades                        = -S_eccentricity:1:S_eccentricity;
    targetOffIntervals{1}           = [stimuliOffset Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(headCenteredTargetLocations);
        
        % Pick saccade
        s = randi(length(saccades));
        
        % Make sure it is big enough and keeps target on retina
        while((abs(saccades(s)) < saccade_threshold) && ((-R_eccentricity + r) <= saccades(s) && saccades(s) <= (r + R_eccentricity)))
            s = randi(length(saccades));
        end
        
        % Location of target
        r = headCenteredTargetLocations(i);
        
        saccadeTargets(i) = saccades(s);
        
        stimuli{i}.headCenteredTargetLocations  = r;
        stimuli{i}.saccadeTargets               = saccades(s);
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.numSaccades                  = length(saccadeOnset);
        
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, 0, saccadeOnset, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;

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
                                    'fixationPeriod', ...
                                    'stimuliOffset', ...
                                    'headCenteredTargetLocations', ...
                                    'saccadeTargets', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end