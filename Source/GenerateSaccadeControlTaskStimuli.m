
%
%  GenerateSaccadeControlTaskStimuli.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateSaccadeControlTaskStimuli(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-SaccadeControlTask'];
    stimulitype = 'SaccadeControlTask';
    
    % Params
    dt                          = 0.010; %(s)
    seed                        = 77;
    S_eccentricity              = 30;
    S_density                   = 1;
    
    % Dynamical quantities
    saccadeOnsetDelay           = 0.1; %(s)
    fixationPeriod              = 0.4; %(s)
    saccadeSpeed                = 300; %(deg/s)

    % Generate stimuli
    rng(seed);
    Duration                    = saccadeOnsetDelay + (2*S_eccentricity/saccadeSpeed) + fixationPeriod; % (s)
    saccadeTargets              = -S_eccentricity:S_density:S_eccentricity;
    targetOffIntervals{1}       = [0 Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(saccadeTargets);
        
        stimuli{i}.initialEyePosition           = 0;
        stimuli{i}.headCenteredTargetLocations  = [];
        stimuli{i}.saccadeTimes                 = saccadeOnsetDelay;
        stimuli{i}.saccadeTargets               = saccadeTargets(i);
        stimuli{i}.numSaccades                  = length(stimuli{i}.saccadeTargets);
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, stimuli{i}.initialEyePosition, stimuli{i}.saccadeTimes, stimuli{i}.saccadeTargets);
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
                                    'S_density', ...
                                    'saccadeTargets', ...
                                    'targetOffIntervals', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'saccadeOnsetDelay', ...
                                    'Duration', ...
                                    'fixationPeriod', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end