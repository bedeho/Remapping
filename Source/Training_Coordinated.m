
%
%  Training_Coordinated.m
%  Remapping
%
%  Created by Bedeho Mender on 04/10/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [Training_RF_Locations, Training_Saccades, filename] = Training_Coordinated(Name, dt)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-Training_Coordinated'];
    stimulitype = 'Training_Coordinated';
    
    % Params
    %dt                              = 0.010; % (s)
    seed                            = 77;
    S_eccentricity                  = 30; % (deg)
    S_density                       = 1; % (deg)
    R_eccentricity                  = 45; % (deg)
    R_density                       = 1; % (deg)
    minimum_Saccade_Amplitude       = 10; % (deg)
    
    % Dynamical quantities
    saccadeSpeed                    = 300; %if changed, then change in GenerateEyeTrace.m as well!
    saccadeOnset                    = 0.200; % w.r.t start of task
    fixationPeriod                  = 0.300; % time from saccade offset

    % Generate stimuli
    rng(seed);
    Duration                        = dtRoundUpPeriod(saccadeOnset + (2*S_eccentricity/saccadeSpeed) + fixationPeriod, dt); % (s), the middle part of sum is to account for maximum saccade times
    saccades                        = -S_eccentricity:S_density:S_eccentricity;
    
    %Training_RF_Locations           = (-R_eccentricity+minimum_Saccade_Amplitude):R_density:(R_eccentricity-minimum_Saccade_Amplitude);
    
    %Training_RF_Locations           = [-20 -17 -15 -12 -10 -7 -5 -2 0 2 5 7 10 12 15 17 20];
    
    %Training_RF_Locations           = [-20 -15 -10 -5 0 5 10 15 20];
    
    %Training_RF_Locations           = [-20 -10 0 10 20];
    
    Training_RF_Locations           = [-20 20];
    
    %Training_RF_Locations           = [0];
    
    hardcoded = false;
    
    %Training_RF_Locations           = [-15 (-15 + 9)]; % (-15 + 9) 
    %Training_Saccades               = [9 -16]; % (-16)
    
    %figure;
    %hold on;

    k = 1;
    for i = 1:length(Training_RF_Locations);
        
        % Receptive field of neuron to be trained
        h = Training_RF_Locations(i);
        
        if(~hardcoded),
        
            % Pick random saccade
            s = randi(length(saccades));

            % Make sure 
            % 1) saccade is big enough
            % 2) FRF is on retina
            % 3) CRF is kept on retina after saccade
            while((abs(saccades(s)) < minimum_Saccade_Amplitude) || ~(-R_eccentricity <= h+saccades(s) && h+saccades(s) <= R_eccentricity) || ~(-R_eccentricity <= h-saccades(s) && h-saccades(s) <= R_eccentricity))
                s = randi(length(saccades));
            end

            % Save saccade
            Training_Saccades(i) = saccades(s);
        end
        
        % FRF Trial
        stimuli{k}.initialEyePosition           = 0;
        stimuli{k}.headCenteredTargetLocations  = h+Training_Saccades(i);
        stimuli{k}.saccadeTimes                 = saccadeOnset;
        stimuli{k}.saccadeTargets               = Training_Saccades(i);
        stimuli{k}.numSaccades                  = length(stimuli{k}.saccadeTargets);
        stimuli{k}.targetOffIntervals           = {[]};

        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
        stimuli{k}.eyePositionTrace             = eyePositionTrace;
        stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
        stimuli{k}.stimOnsetTimes               = 0;
        
        %plot(stimuli{k}.headCenteredTargetLocations, stimuli{k}.saccadeTargets, 'or');

        k = k + 1;
        
        %{
        % CRF Trial
        stimuli{k}.initialEyePosition           = 0;
        stimuli{k}.headCenteredTargetLocations  = h;%-Training_Saccades(i);
        stimuli{k}.saccadeTimes                 = saccadeOnset;
        stimuli{k}.saccadeTargets               = Training_Saccades(i);
        stimuli{k}.numSaccades                  = length(stimuli{k}.saccadeTargets);
        stimuli{k}.targetOffIntervals           = {[]};

        [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
        stimuli{k}.eyePositionTrace             = eyePositionTrace;
        stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
        stimuli{k}.stimOnsetTimes               = 0;
        
        plot(stimuli{k}.headCenteredTargetLocations, stimuli{k}.saccadeTargets, 'ob');
        
        k = k + 1;
        %}
        
    end
    
    
    plot(Training_RF_Locations, Training_Saccades, 'or');
    xlabel('RF Location (deg)');
    xlim([-R_eccentricity R_eccentricity]);
    ylabel('Saccade (deg)');
    ylim([-S_eccentricity S_eccentricity]);
    
   
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
                                    'saccadeOnset', ...
                                    'fixationPeriod', ...
                                    'stimulitype', ...
                                    'stimuli', ...
                                    'Duration', ...
                                    'saccadeSpeed', ...
                                    'dt', ...
                                    'seed');
                                
                                

end