
%
%  Testing_SaccadeControl.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_SaccadeControl(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-SaccadeControl'];
    stimulitype = 'SaccadeControl';
    
    % Params
    dt                          = 0.010; %(s)
    seed                        = 77;
    S_eccentricity              = 30;
    S_density                   = 1; % should be 1, we test everyones response
    
    % Dynamical quantities
    saccadeOnset                = 0.1; %(s) w.r.t start of task
    postSaccadefixationPeriod   = 0.4; %(s)
    saccadeSpeed                = 300; %(deg/s)

    % Generate stimuli
    rng(seed);
    Duration                    = saccadeOnset + postSaccadefixationPeriod; % (s)
    saccadeTargets              = -S_eccentricity:S_density:S_eccentricity;
    targetOffIntervals{1}       = [0 Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(saccadeTargets);
        
        stimuli{i}.initialEyePosition           = 0;
        stimuli{i}.headCenteredTargetLocations  = [];
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.saccadeTargets               = saccadeTargets(i);
        stimuli{i}.numSaccades                  = length(stimuli{i}.saccadeTargets);
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, stimuli{i}.initialEyePosition, stimuli{i}.saccadeTimes, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;
        stimuli{i}.stimOnsetTimes               = [];
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
                                    'saccadeTargets', ...
                                    'targetOffIntervals', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'saccadeOnset', ...
                                    'Duration', ...
                                    'postSaccadefixationPeriod', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end