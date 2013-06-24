
%
%  GenerateStimuliControlTaskStimuli.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateStimuliControlTaskStimuli(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-StimuliControlTask'];
    stimulitype = 'StimuliControlTask';
    
    % Params
    dt                          = 0.010; %(s)
    seed                        = 77;
    R_eccentricity              = 45;
    R_density                   = 1;
    
    % Dynamical quantities
    stimuliOnsetDelay           = 0.1; %(s)
    stimuliDuration             = 0.1; %(s)
    stimuliOffsetPeriod         = 0.3; %(s)

    % Generate stimuli
    rng(seed);
    Duration                    = stimuliOnsetDelay+stimuliDuration+stimuliOffsetPeriod; % (s)
    headCenteredTargetLocations = -R_eccentricity:R_density:R_eccentricity;
    targetOffIntervals{1}       = [0 stimuliOnsetDelay;(stimuliOnsetDelay+stimuliDuration) Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(headCenteredTargetLocations);
        
        stimuli{i}.initialEyePosition           = 0;
        stimuli{i}.headCenteredTargetLocations  = headCenteredTargetLocations(i);
        stimuli{i}.saccadeTimes                 = [];
        stimuli{i}.saccadeTargets               = [];
        stimuli{i}.numSaccades                  = length(stimuli{i}.saccadeTargets);
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        stimuli{i}.eyePositionTrace             = GenerateEyeTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, stimuli{i}.initialEyePosition, stimuli{i}.saccadeTimes, stimuli{i}.saccadeTargets);
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    mkdir(stimuliFolder);
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'stimulitype', ...
                                    'headCenteredTargetLocations', ...
                                    'targetOffIntervals', ...
                                    'R_eccentricity', ...
                                    'R_density', ...
                                    'Duration', ...
                                    'stimuli', ...
                                    'stimuliOnsetDelay', ...
                                    'stimuliDuration', ...
                                    'stimuliOffsetPeriod', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end