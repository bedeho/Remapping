
%
%  Training_LHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 13/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function [Training_RF_Locations, Training_Saccades, filename] = Training_LHeiser(Name, dt)

    % Import global variables
    declareGlobalVars();
    global base;

    filename = [Name '-Training_LHeiser'];
    stimulitype = 'Training_LHeiser';
    
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
    
    %Unique_Training_RF_Locations           = (-R_eccentricity+minimum_Saccade_Amplitude):R_density:(R_eccentricity-minimum_Saccade_Amplitude);
    
    %Unique_Training_RF_Locations           = [-20 -17 -15 -12 -10 -7 -5 -2 0 2 5 7 10 12 15 17 20];
    
    Unique_Training_RF_Locations           = [-20 20];
    
    numberOfDirections              = 3;

    k = 1;
    for i = 1:length(Unique_Training_RF_Locations);
        
        % Receptive field of neuron to be trained
        h = Unique_Training_RF_Locations(i);
        
        for d = 1:numberOfDirections,
        
            % Pick random saccade
            s = randi(length(saccades));

            % Make sure 
            % 1) saccade is big enough
            % 2) FRF is on retina
            % 3) CRF is kept on retina after saccade
            while((abs(saccades(s)) < minimum_Saccade_Amplitude) || ~(-R_eccentricity <= h+saccades(s) && h+saccades(s) <= R_eccentricity) || ~(-R_eccentricity <= h-saccades(s) && h-saccades(s) <= R_eccentricity))
                s = randi(length(saccades));
            end

            % Save saccade and rf
            Training_Saccades(k)        = saccades(s);
            Training_RF_Locations(k)    = h;

            % FRF Trial
            stimuli{k}.initialEyePosition           = 0;
            stimuli{k}.headCenteredTargetLocations  = h+Training_Saccades(k);
            stimuli{k}.saccadeTimes                 = saccadeOnset;
            stimuli{k}.saccadeTargets               = Training_Saccades(k);
            
            stimuli{k}.numSaccades                  = length(stimuli{k}.saccadeTargets);
            stimuli{k}.targetOffIntervals           = {[]};

            [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, stimuli{k}.headCenteredTargetLocations, stimuli{k}.targetOffIntervals, stimuli{k}.initialEyePosition, stimuli{k}.saccadeTimes, stimuli{k}.saccadeTargets);
            stimuli{k}.eyePositionTrace             = eyePositionTrace;
            stimuli{k}.retinalTargetTraces          = retinalTargetTraces;
            stimuli{k}.stimOnsetTimes               = 0;


            k = k + 1;
        
        end

        
    end
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep filename];
    
    if exist(stimuliFolder),
        system(['rm -R ' stimuliFolder]);
    end 
    
    mkdir(stimuliFolder);
    
    % Save stimuli
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
                                    'numberOfDirections', ...
                                    'dt', ...
                                    'seed');
                                
    % Make, save & close figure
    f = figure('Units','pixels','position', [1000 1000 420 300]);
    plot(Training_RF_Locations, Training_Saccades, 'or');

    hXLabel = xlabel('Retinal Locaton (deg)');
    hYLabel = ylabel('Saccade (deg)');
    set([hYLabel hXLabel], 'FontSize', 14);
    set(gca, 'FontSize', 12);
    xlim(1.1*[-R_eccentricity R_eccentricity]);
    ylim(1.1*[-S_eccentricity S_eccentricity]);
    pbaspect([2*S_eccentricity 2*S_eccentricity 1]);

    saveas(f, [stimuliFolder filesep 'trainingstimuli.eps']);
    close(f);

end