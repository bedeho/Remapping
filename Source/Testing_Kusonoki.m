
%
%  Testing_Kusonoki.m
%  Remapping
%
%  Created by Bedeho Mender on 23/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Testing_Kusonoki(Name, dt)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-Kusonoki'];
    stimulitype = 'Kusonoki';
    
    %% Parameters
    
    % Technical
    %dt                          = 0.010;
    seed                            = 77;
    rng(seed);
    
    % Spatial
    S_eccentricity                  = 30;
    R_eccentricity                  = 45;
    saccade_threshold               = S_eccentricity/2;
    
    % Temporal
    saccadeSpeed                    = 300; % (deg/s) if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.600; % w.r.t start of task    
    stimulusDuration                = 0.100;
    stimulusOnsetTimes              = saccadeOnset + ((-(0.400+stimulusDuration)):0.050:0.100); % w.r.t start of trial

    % Utilities - derived
    screen_locations                = [-10];% -5 0 5 10];%-R_eccentricity:1:R_eccentricity; %-7
    saccades                        = -S_eccentricity:1:S_eccentricity;
    saccadeDelayTime                = roundn((2*S_eccentricity/saccadeSpeed)+0.05,-1); % round to nearest hundred above
    Duration                        = max(stimulusOnsetTimes(end),saccadeOnset + saccadeDelayTime) + 0.400; % (s), make sure we have enough time after last stim onset time to have
    % space for response window!
    
    % z=1: Current RF trials
    % z=2: Future RF trials

    % ITerate location where rf will be, current or future.
    trialNr = 1;
    for i = 1:length(screen_locations),

        % Get rf location
        RF_location = screen_locations(i);

        % Pick a first random saccade
        s = randi(length(saccades));

        % Make sure this saccade is big enough, and that future RF is
        % on retina
        while((abs(saccades(s)) < saccade_threshold) || ~(-R_eccentricity <= RF_location+saccades(s) && RF_location+saccades(s) <=R_eccentricity))
            s = randi(length(saccades));
        end

        for z=1:2,

            % Where to put stimuli prior to saccade
            if(z==1), %current rf trial
                stim_location = RF_location;
            else
                stim_location = RF_location + saccades(s);
            end

            % Iterate different stimulus onset times
            for t = 1:length(stimulusOnsetTimes),

                onsetTime                               = stimulusOnsetTimes(t);
                offsetTime                              = onsetTime + stimulusDuration;
                targetOffIntervals{1}                   = [0 (onsetTime-dt);offsetTime Duration]; % (s) [start_OFF end_OFF; start_OFF end_OFF]

                % Generate trace
                stimuli{trialNr}.headCenteredTargetLocations  = stim_location;
                stimuli{trialNr}.saccadeTargets               = saccades(s);
                stimuli{trialNr}.saccadeTimes                 = saccadeOnset;
                stimuli{trialNr}.targetOffIntervals           = targetOffIntervals;
                [eyePositionTrace, retinalTargetTraces]       = GenerateTrace(Duration, dt, stimuli{trialNr}.headCenteredTargetLocations, stimuli{trialNr}.targetOffIntervals, 0, stimuli{trialNr}.saccadeTimes, stimuli{trialNr}.saccadeTargets);
                stimuli{trialNr}.eyePositionTrace             = eyePositionTrace;
                stimuli{trialNr}.retinalTargetTraces          = retinalTargetTraces;
                stimuli{trialNr}.stimOnsetTimes               = onsetTime;

                % Add simple information
                stimuli{trialNr}.neuron_RF_location           = RF_location;
                stimuli{trialNr}.trialType                    = z;
                stimuli{trialNr}.targetNr                     = i;
                stimuli{trialNr}.saccadeNr                    = s;
                stimuli{trialNr}.stimOnsetNr                  = t;

                trialNr = trialNr + 1;
            end

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
                                    'R_eccentricity', ...
                                    'saccadeOnset', ...
                                    'stimulusDuration', ...
                                    'stimulusOnsetTimes', ...
                                    'saccades', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end