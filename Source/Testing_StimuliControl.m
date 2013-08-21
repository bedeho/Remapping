
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
    R_density                   = 1; % should be 1, we test everyones latency
    
    % Dynamical quantities
    stimuliOnset           = 0.100; %(s)
    stimuliDuration             = 0.300; %(s)
    stimuliOffsetPeriod         = 0.300; %(s)

    % Generate stimuli
    rng(seed);
    Duration                    = stimuliOnset+stimuliDuration+stimuliOffsetPeriod; % (s)
    headCenteredTargetLocation  = -R_eccentricity:R_density:R_eccentricity;
    targetOffIntervals{1}       = [0 stimuliOnset;(stimuliOnset+stimuliDuration) Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]
    
    for i = 1:length(headCenteredTargetLocation);
        
        stimuli{i}.initialEyePosition           = 0;
        stimuli{i}.headCenteredTargetLocations   = headCenteredTargetLocation(i);
        stimuli{i}.saccadeTimes                 = [];
        stimuli{i}.saccadeTargets               = [];
        stimuli{i}.numSaccades                  = length(stimuli{i}.saccadeTargets);
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, targetOffIntervals, stimuli{i}.initialEyePosition, stimuli{i}.saccadeTimes, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;
        stimuli{i}.stimOnsetTimes               = stimuliOnset;
        
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'stimulitype', ...
                                    'targetOffIntervals', ...
                                    'R_eccentricity', ...
                                    'R_density', ...
                                    'Duration', ...
                                    'stimuli', ...
                                    'stimuliOnset', ...
                                    'stimuliDuration', ...
                                    'stimuliOffsetPeriod', ...
                                    'dt', ...
                                    'seed');
end