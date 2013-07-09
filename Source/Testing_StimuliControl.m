
%
%  Testing_StimuliControl.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_StimuliControl(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-StimuliControl'];
    stimulitype = 'StimuliControl';
    
    % Params
    dt                          = 0.010; %(s)
    seed                        = 77;
    R_eccentricity              = 45;
    R_density                   = 10;
    
    % Dynamical quantities
    stimuliOnsetDelay           = 0.100; %(s)
    stimuliDuration             = 0.300; %(s)
    stimuliOffsetPeriod         = 0.300; %(s)

    % Generate stimuli
    rng(seed);
    Duration                    = stimuliOnsetDelay+stimuliDuration+stimuliOffsetPeriod; % (s)
    headCenteredTargetLocations = [0];%-R_eccentricity:R_density:R_eccentricity;
    targetOffIntervals{1}       = [0 stimuliOnsetDelay;(stimuliOnsetDelay+stimuliDuration) Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(headCenteredTargetLocations);
        
        stimuli{i}.initialEyePosition           = 0;
        stimuli{i}.headCenteredTargetLocations  = headCenteredTargetLocations(i);
        stimuli{i}.saccadeTimes                 = [];
        stimuli{i}.saccadeTargets               = [];
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
                                    'dt', ...
                                    'seed');
end