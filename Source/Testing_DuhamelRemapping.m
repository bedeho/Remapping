
%
%  Testing_DuhamelRemapping.m
%  Remapping
%
%  Created by Bedeho Mender on 06/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_DuhamelRemapping(Name, dt, Training_RF_Locations, Training_Saccades, stimulitype, saccadeOnset, stimuliDuration, stimuliOnset, postSaccadefixationPeriod)

    % Import global variables
    declareGlobalVars();
    global base;
    
    %% Parameters
    
    % Technical
    %dt                          = 0.010;
    seed                            = 77;
    rng(seed);
    
    % Spatial
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    saccade_threshold               = S_eccentricity/2;
    
    % Dynamical
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    
    % Training_RF_Locations, Training_Saccades
    
    if nargin < 9,
        
        % Temporal
        stimulitype                     = 'DuhamelRemapping';
        saccadeOnset                    = 0.300; % (s) w.r.t start of task
        stimuliDuration                 = 0.100; % 0.100; (s)
        stimuliOnset                    = saccadeOnset - 2*stimuliDuration; % (s) w.r.t start of task
        postSaccadefixationPeriod       = 0.300; % (s) time from saccade COMPLETION ESTIMATE "saccadeDelayTime" below.
            
        if nargin < 4,
            
            % Utilities - derived
            currentRF = 20%0%8; %[-10 -5 0 5 10];%-R_eccentricity:1:R_eccentricity; % 10; Remapping TARGET, i.e. postsaccadic (-R_eccentricity+R_edge_effect_buffer):1:(R_eccentricity+R_edge_effect_buffer)
            saccades  = -20%-S_eccentricity:1:S_eccentricity; % Pick among these saccades
        else
            currentRF = Training_RF_Locations;
        end

    else
        currentRF = Training_RF_Locations;
    end
    
    saccadeDelayTime                = (2*S_eccentricity/saccadeSpeed) + 0.05; % round to nearest hundred above
    Duration                        = dtRoundUpPeriod(saccadeOnset + saccadeDelayTime + postSaccadefixationPeriod, dt); % (s), the middle part of sum is to account for maximum saccade times
    
    % classic
    targetOffIntervals{1}           = [0 (stimuliOnset-dt); (stimuliOnset+stimuliDuration) Duration];
    %%targetOffIntervals{1}           = [0 stimuliOnset; (stimuliOnset+stimuliDuration) Duration];
    
    %% Generate stimuli
    for i = 1:length(currentRF);
        
        % Rf that will get remapping activity INTO it
        r = currentRF(i);
        
        if nargin < 4,
            
            % Pick saccade
            s = randi(length(saccades));

            % Make sure this saccade is big enough, and has current_ref on retina
            while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= r+saccades(s) && r+saccades(s) <=R_eccentricity))
                s = randi(length(saccades));
            end

            saccade = saccades(s);
        else
            saccade = Training_Saccades(i);
        end
        
        % Screen location of stimulus
        futureRF = r+saccade;
        
        % Generate trace
        stimuli{i}.headCenteredTargetLocations  = futureRF;
        stimuli{i}.saccadeTargets               = saccade;
        stimuli{i}.saccadeTimes                 = saccadeOnset;
        stimuli{i}.targetOffIntervals           = targetOffIntervals;
        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{i}.headCenteredTargetLocations, stimuli{i}.targetOffIntervals, 0, saccadeOnset, stimuli{i}.saccadeTargets);
        stimuli{i}.eyePositionTrace             = eyePositionTrace;
        stimuli{i}.retinalTargetTraces          = retinalTargetTraces;
        stimuli{i}.stimOnsetTimes               = stimuliOnset;
        
        % Add simple information
        stimuli{i}.futureRF  = futureRF;
        stimuli{i}.currentRF = r;
    end
    
    % Save params
    filename = [Name '-' stimulitype];
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'S_eccentricity', ...
                                    'R_eccentricity', ...
                                    'saccadeOnset', ...
                                    'stimuliDuration', ...
                                    'stimuliOnset', ...
                                    'postSaccadefixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end