
%
%  Testing_DuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 06/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_DuhamelRemapping(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-DuhamelRemapping'];
    stimulitype = 'DuhamelRemapping';
    
    % Params
    dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    
    saccade_threshold               = S_eccentricity/4;
    
    % Dynamical quantities
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.100; % (s) w.r.t start of task
    fixationPeriod                  = 0.300; % (s) time from saccade onset
    
    % Generate stimuli
    rng(seed);
    Duration                        = saccadeOnset + fixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    headCenteredTargetLocations     = -R_eccentricity:1:R_eccentricity;
    saccades                        = -S_eccentricity:1:S_eccentricity;
    
    for i = 1:length(headCenteredTargetLocations);
        
        % Location of target
        r = headCenteredTargetLocations(i);
        
        % Pick saccade
        s = randi(length(saccades));
        
        while((abs(saccades(s)) > saccade_threshold)) %(-R_eccentricity + r <= saccades(s)) && (saccades(s) <= r + R_eccentricity)
            s = randi(length(saccades));
        end
        
        stimuli{i}.headCenteredTargetLocations  = r + saccades(s);
        stimuli{i}.saccadeTarget                = saccades(s);
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, [], 0, saccadeOnset, stimuli{k}.saccadeTargets);
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
                                    'headCenteredTargetLocations', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end