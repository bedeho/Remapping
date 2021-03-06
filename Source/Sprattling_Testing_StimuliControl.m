
%
%  Sprattling_Testing_StimuliControl.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Sprattling_Testing_StimuliControl(Name, dt) %, headCenteredTargetLocation)

    % Import global variables
    declareGlobalVars();
    global base;

    % Params
    %dt                          = 0.010; %(s)
    seed                        = 77;
    R_eccentricity              = 45;
    R_density                   = 1; % should be 1, we test everyones latency
    
    filename = [Name '-Sprattling_StimuliControl'];
    stimulitype = 'Sprattling_StimuliControl';
    headCenteredTargetLocation  = -38; %;-R_eccentricity:R_density:R_eccentricity;
        
    %{
    if nargin < 3,
        filename = [Name '-StimuliControl'];
        stimulitype = 'StimuliControl';
        headCenteredTargetLocation  = -R_eccentricity:R_density:R_eccentricity;
    else
        filename = [Name '-StimuliControl'];
        stimulitype = 'StimuliControl2';
    end
    %}
    
    % Dynamical quantities
    stimuliOnset                = 0.100; %(s) classic = 0.100
    stimuliDuration             = 0.100; %(s) classic = 0.300
    stimuliOffsetPeriod         = 0.400+0.300; %(s) classic = 0.400, LHeiser2005 = offset before saccade + sacade onset aligned 300 ms window = 0.400+0.300

    % Generate stimuli
    rng(seed);
    Duration                    = dtRoundUpPeriod(stimuliOnset+stimuliDuration+stimuliOffsetPeriod, dt); % (s)
    
    targetOffIntervals{1}       = [0 (stimuliOnset-dt);(stimuliOnset+stimuliDuration) Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]

    
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