
%
%  Testing_Kusonoki.m
%  Remapping
%
%  Created by Bedeho Mender on 23/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_Kusonoki(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-Kusonoki'];
    stimulitype = 'Kusonoki';
    
    % Params
    dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30;
    S_density                       = 10;
    R_eccentricity                  = 45;
    R_density                       = 10;
    
    % Dynamical quantities
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.4; % w.r.t start of task
    earliestStimulusOnsetTime       = 0.100; % w.r.t start of task
    lastStimulusOnsetTime           = saccadeOnset + 0.500; % w.r.t start of task
    stimulusDuration                = 0.1;
    fixationPeriod                  = lastStimulusOnsetTime + stimulusDuration + 0.100; % time from saccade onset
    stimulusOnsetTimes              = earliestStimulusOnsetTime:0.200:lastStimulusOnsetTime;% earliestStimulusOnsetTime

    % Generate stimuli
    rng(seed);
    Duration                        = saccadeOnset + fixationPeriod; % (s), the middle part of sum is to account for maximum saccade times
    saccadeTargets                  = -S_eccentricity:S_density:S_eccentricity;
    headCenteredTargetLocations     = -R_eccentricity:R_density:R_eccentricity;

    k = 1;
    for i = 1:length(headCenteredTargetLocations);
        
        % Location of target
        h = headCenteredTargetLocations(i);
        
        % Deduce all valid saccades keeping this target on retina
        sTargets = saccadeTargets(-R_eccentricity + h <= saccadeTargets <= h + R_eccentricity);
        
        for j = 1:length(sTargets),
            
            % Get saccade
            s = sTargets(j);
            
            for t = 1:length(stimulusOnsetTimes),
                
                    onsetTime                               = stimulusOnsetTimes(t);
                    offsetTime                              = onsetTime + stimulusDuration;
                    targetOffIntervals{1}                   = [0 onsetTime;offsetTime Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
            
                    stimuli{k}.initialEyePosition           = 0;
                    stimuli{k}.headCenteredTargetLocations  = h;
                    stimuli{k}.saccadeTimes                 = saccadeOnset;
                    stimuli{k}.saccadeTargets               = s;
                    stimuli{k}.numSaccades                  = length(stimuli{k}.saccadeTargets);
                    stimuli{k}.targetOffIntervals           = targetOffIntervals;
                    
                    
                    [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
                    stimuli{k}.eyePositionTrace             = eyePositionTrace;
                    stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
                    
                    % meta data
                    stimuli{k}.targetNr                     = i;
                    stimuli{k}.saccadeNr                    = j;
                    stimuli{k}.stimOnsetNr                  = t;

                    k = k + 1;
            end
            
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
                                    'earliestStimulusOnsetTime', ...
                                    'lastStimulusOnsetTime', ...
                                    'stimulusDuration', ...
                                    'fixationPeriod', ...
                                    'stimulusOnsetTimes', ...
                                    'saccadeTargets', ...
                                    'headCenteredTargetLocations', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end