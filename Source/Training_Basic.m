
%
%  Training_Basic.m
%  Remapping
%
%  Created by Bedeho Mender on 23/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Training_Basic(Name,dt)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-Training'];
    stimulitype = 'Training';
    
    % Params
    %dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30;
    S_density                       = 5;
    R_eccentricity                  = 45;
    R_density                       = 5;
    
    % Dynamical quantities
    saccadeSpeed                    = 300; %if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.200; % w.r.t start of task
    fixationPeriod                  = 0.300; % time from saccade offset

    % Generate stimuli
    rng(seed);
    Duration                        = saccadeOnset + (2*S_eccentricity/saccadeSpeed) + fixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    saccadeTargets                  = -S_eccentricity:S_density:S_eccentricity; %[-20]%
    headCenteredTargetLocations     = -R_eccentricity:R_density:R_eccentricity;

    k = 1;
    for i = 1:length(headCenteredTargetLocations);
        
        % Location of target
        h = headCenteredTargetLocations(i);
        
        % Deduce all valid saccades keeping this target on retina
        sTargets = saccadeTargets((-R_eccentricity + h <= saccadeTargets) & (saccadeTargets <= h + R_eccentricity));
        
        for j = 1:length(sTargets),
            
            % Get saccade
            s = sTargets(j);
            
            stimuli{k}.initialEyePosition           = 0;
            stimuli{k}.headCenteredTargetLocations  = h;
            stimuli{k}.saccadeTimes                 = saccadeOnset;
            stimuli{k}.saccadeTargets               = s;
            stimuli{k}.numSaccades                  = length(stimuli{k}.saccadeTargets);
            stimuli{k}.targetOffIntervals           = {[]};

            [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
            stimuli{k}.eyePositionTrace             = eyePositionTrace;
            stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
            stimuli{k}.stimOnsetTimes               = 0;
            
            k = k + 1;
            
        end
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'S_eccentricity', ...
                                    'S_density', ...
                                    'R_eccentricity', ...
                                    'R_density', ...
                                    'saccadeOnset', ...
                                    'fixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end