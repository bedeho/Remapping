
%
%  Testing_CLayerProbeTask.m
%  Remapping
%
%  Created by Bedeho Mender on 13/08/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_CLayerProbeTask(Name, dt)

    % Import global variables
    declareGlobalVars();
    global base;
    
    if(nargin < 1)
        Name = 'basic';
    end

    filename = [Name '-CLayerProbe'];
    stimulitype = 'CLayerProbe';
    
    % Params
    %dt                          = 0.010;
    seed                        = 77; 
    S_eccentricity              = 30;
    S_density                   = 4;
    R_eccentricity              = 45;
    R_density                   = 5;

    % Generate stimuli
    rng(seed);
    saccadeOnset                    = 0.200; % (s)
    Duration                        = dtRoundUpPeriod(saccadeOnset+0.100, dt); % (s), recording window is 50ms after saccade, so atleast this much delay is needed,the middle part of sum is to account for maximum saccade times
    saccadeTargets                  = -S_eccentricity:S_density:S_eccentricity;
    headCenteredTargetLocations     = -R_eccentricity:R_density:R_eccentricity;

    k = 1;
    for i = 1:length(headCenteredTargetLocations),
        
        % Location of target
        h = headCenteredTargetLocations(i);
            
        for j = 1:length(saccadeTargets),
            
            % Saccade target
            s = saccadeTargets(j);

            stimuli{k}.initialEyePosition           = 0;
            stimuli{k}.headCenteredTargetLocations  = h;
            stimuli{k}.saccadeTimes                 = saccadeOnset;
            stimuli{k}.saccadeTargets               = s;
            stimuli{k}.numSaccades                  = 1;
            stimuli{k}.targetOffIntervals           = {[]};

            [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
            stimuli{k}.eyePositionTrace             = eyePositionTrace;
            stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
            stimuli{k}.stimOnsetTimes               = 0;

            % meta data
            stimuli{k}.targetNr                     = i;
            stimuli{k}.saccadeNr                    = j;

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
                                    'saccadeTargets', ...
                                    'headCenteredTargetLocations', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'saccadeOnset', ...
                                    'Duration', ...
                                    'dt', ...
                                    'seed');
end