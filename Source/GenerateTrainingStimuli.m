
%
%  GenerateTrainingStimuli.m
%  Remapping
%
%  Created by Bedeho Mender on 23/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateTrainingStimuli(Name)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-GenerateTrainingStimuli'];
    stimulitype = 'TrainingStimuli';
    
    % Params
    dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30;
    S_density                       = 1;
    R_eccentricity                  = 45;
    R_density                       = 1;
    
    % Dynamical quantities
    saccadeOnset                    = 0.4; % w.r.t start of task
    earliestStimulusOnsetTime       = 0.100; % w.r.t start of task
    lastStimulusOnsetTime           = saccadeOnset + 0.500; % w.r.t start of task
    stimulusDuration                = 0.1;
    fixationPeriod                  = lastStimulusOnsetTime + stimulusDuration + 100; % time from saccade onset
    stimulusOnsetTimes              = earliestStimulusOnsetTime:100:lastStimulusOnsetTime;
    
    % Saccadic
    saccadeSpeed                    = 300; % (deg/s)

    % Generate visual target locations
    rng(seed);
    Duration                        = saccadeOnsetDelay + (saccadeSpeed/2*S_eccentricity) + fixationPeriod; % (s)
    saccadeTargets                  = -S_eccentricity:S_density:S_eccentricity;
    headCenteredTargetLocations     = -R_eccentricity:R_density:R_eccentricity;
    
    k = 1;
    for i = 1:length(headCenteredTargetLocations);
        
        % Location of target
        h = headCenteredTargetLocations(i);
        
        % Deduce all valid saccades keeping this target on retina
        saccadeTargets = saccadeTargets(-R_eccentricity + h <= saccadeTargets <= h + R_eccentricity);
        
        for j = 1:length(saccadeTargets),
            
            % Get saccade
            s = saccadeTargets(j);
            
            for t = 1:length(stimulusOnsetTimes),
                
                    onsetTime = stimulusOnsetTimes(t);
            
                    stimuli{k}.initialEyePosition = 0;
                    stimuli{k}.headCenteredTargetLocations = h;
                    stimuli{k}.saccadeTimes = saccadeOnset;
                    stimuli{k}.saccadeTargets = s;
                    stimuli{k}.numSaccades = length(stimuli{k}.saccadeTargets);

                    stimuli{k}.eyePositionTrace = GenerateEyeTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, {[]}, 0, saccadeSpeed, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);

                    k = k + 1;
            end
            
        end
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    mkdir(stimuliFolder);
    save([stimuliFolder filesep 'stim.mat'] , ...
                                    'S_eccentricity', ...
                                    'S_density', ...
                                    'R_eccentricity', ...
                                    'R_density', ...
                                    'saccadeOnset', ...
                                    'earliestStimulusOnsetTime', ...
                                    'lastStimulusOnsetTime', ...
                                    'stimulusDuration', ...
                                    'fixationPeriod', ...
                                    'stimulusOnsetTimes', ...
                                    'saccadeTargets', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
end