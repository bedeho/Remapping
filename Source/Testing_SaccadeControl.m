
%
%  Testing_SaccadeControl.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_SaccadeControl(Name, dt, saccadeTargets)

    % Import global variables
    declareGlobalVars();
    global base;
    
    % Params
    %dt                          = 0.010;
    seed                        = 77;
    S_eccentricity              = 30;
    S_density                   = 1; % should be 1, we test everyones response
    
    filename       = [Name '-SaccadeControl'];
    stimulitype    = 'SaccadeControl';
    saccadeTargets = -S_eccentricity:S_density:S_eccentricity;
    
    %{    
    if nargin<3,
        filename       = [Name '-SaccadeControl'];
        stimulitype    = 'SaccadeControl';
        saccadeTargets = -S_eccentricity:S_density:S_eccentricity;
    else
        filename       = [Name '-SaccadeControl'];
        stimulitype    = 'SaccadeControl2';
    end
    %}
        
    % Dynamical quantities
    saccadeOnset                = 0.100; %(s) w.r.t start of task, CLASSIC 0.150, to make identical to stim control
    postSaccadeOnsetPeriod      = 0.800; % CLASSIC: postSaccadefixationPeriod   = 0.400; %(s)
    %saccadeSpeed                = 300; %(deg/s)

    % Generate stimuli
    rng(seed);
    Duration                    = dtRoundUpPeriod(saccadeOnset + postSaccadeOnsetPeriod, dt); % dtRoundUpPeriod(saccadeOnset + postSaccadefixationPeriod, dt); % (s)
    
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
                                    'postSaccadeOnsetPeriod', ...
                                    'dt', ...
                                    'seed');
end