
%
%  Testing_LHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 13/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function Testing_LHeiser(Name, dt, Training_RF_Locations)

    % Import global variables
    declareGlobalVars();
    global base;
    
    %% Parameters
    
    % Technical
    seed                            = 77;
    rng(seed);
    
    % Spatial
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    saccade_threshold               = S_eccentricity/2;
    
    % Dynamical
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    
    % Training_RF_Locations, Training_Saccades
    saccadeDelayTime                = (2*S_eccentricity/saccadeSpeed) + 0.050; % round to nearest hundred above
     
    % Temporal
    stimulitype                     = 'LHeiser';
    saccadeOnset                    = 0.300; % (s) w.r.t start of task
    stimuliOnset                    = 0.100; % saccadeOnset - 2*stimuliDuration; % (s) w.r.t start of task
    postSaccadefixationPeriod       = 0.300; % (s) time from saccade COMPLETION ESTIMATE "saccadeDelayTime" below.
    stimuliDuration                 = 0.050; % (s)
    
    numberOfDirections              = 3;
    
    currentRF                       = Training_RF_Locations;
    Duration                        = dtRoundUpPeriod(saccadeOnset + saccadeDelayTime + postSaccadefixationPeriod, dt); % (s), the middle part of sum is to account for maximum saccade times

    % classic
    targetOffIntervals{1}           = [0 (stimuliOnset-dt); (stimuliOnset+stimuliDuration) Duration];
    
    %% Generate stimuli
    k = 1;
    for i = 1:length(currentRF);
        
        % Rf that will get remapping activity INTO it
        r = currentRF(i);
        
        for d = 1:numberOfDirections,
        
            % Pick saccade
            s = randi(length(saccades));

            % Make sure this saccade is big enough, and has current_ref on retina
            while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= r+saccades(s) && r+saccades(s) <=R_eccentricity))
                s = randi(length(saccades));
            end

            saccade = saccades(s);

            % Screen location of stimulus
            futureRF = r+saccade;

            % Generate trace
            stimuli{k}.headCenteredTargetLocations  = futureRF;
            stimuli{k}.saccadeTargets               = saccade;
            stimuli{k}.saccadeTimes                 = saccadeOnset;
            stimuli{k}.targetOffIntervals           = targetOffIntervals;
            [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, 0, saccadeOnset, stimuli{k}.saccadeTargets);
            stimuli{k}.eyePositionTrace             = eyePositionTrace;
            stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
            stimuli{k}.stimOnsetTimes               = stimuliOnset;
            stimuli{k}.directionNr                  = d;
            
            % Add simple information
            stimuli{k}.futureRF  = futureRF;
            stimuli{k}.currentRF = r;
            
            % Next 
            k = k + 1;
        
        end
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