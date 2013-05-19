
%
%  GenerateTrainingStimuli.m
%  Remapping
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateTrainingStimuli()

    % Import global variables
    declareGlobalVars();
    global base;

    Name = 'trainingstimuli';

    % Dynamical quantities
    Duration = 2; % (s)
    dt = 0.010; % (s)
    numTimeSteps = ceil(Duration/dt);

    % Random seed
    seed = 77;
    rng(seed);

    % Visual
    headCenteredTargetLocations = [15]; % (deg)
    maxNumberOfVisibleTargets = length(headCenteredTargetLocations);
    
    targetOffIntervals{1} = [];% [0.6 1.8]; %[0.5 0.9;]; % (s) [start_OFF end_OFF; start_OFF end_OFF] <==== Make dt multiples
    %targetOffIntervals{2} = []; %[0.1 0.2;];
    
    assert(length(targetOffIntervals) >= maxNumberOfVisibleTargets, 'On off history not provided for all targets.');
    
    % Saccadic
    initialEyePosition = 0; % (deg), init = 13 deg
    saccadeSpeed = 300; % (deg/s)
    saccadeTimes = [1.2]; % (s) % <==== Make dt multiples, ALLOW FOR SUFFIEIENCT INTERSACCADE TIME TO COMPLETE SACCADES!
    saccadeTargets = [20]; % (deg)
    numSaccades = length(saccadeTimes);
    
    assert(length(saccadeTimes) >= length(saccadeTargets), 'Number of saccade times and targets must match.');

    % Allocate space for traces
    % Index i => Time (i-1)*dt
    eyePositionTrace = zeros(1, numTimeSteps);
    eyePositionTrace(1) = initialEyePosition;

    % Generate eye position trace
    numCompletedSaccades = 0;
    for t=2:numTimeSteps,
        
        presentTime = stepToTime(t);
        
        % Have we completed all saccades beginning before time timestep t?
        if numCompletedSaccades == numSaccades || presentTime <= saccadeTimes(numCompletedSaccades+1),

            % yes, so lets just continue fixating
            eyePositionTrace(t) = eyePositionTrace(t-1);

        else
            % no, there was a not completed saccade beginning before
            % timestep t
            
            % was this the first step across this saccade onset time?
            previousTimeSteptime = stepToTime(t-1);
            timeOffset = abs(previousTimeSteptime - saccadeTimes(numCompletedSaccades+1));
            
            if(timeOffset < dt),
                
                % yes it was, so mini saccade starting
                % point needs to correct for delay between prior time step
                % and saccade onset event to find saccade start position
                
                reduceSaccadeTimeInPresentTimeStepWith = timeOffset;
            else
                
                % no, we passed at some previous time, so mini saccade starting
                % point is just last eye position
                
                reduceSaccadeTimeInPresentTimeStepWith = 0;
            end
            
            % Will one more stime step saccading from
            % miniSaccadeStartPosition take us past saccade target?
            
            saccadeOffset = eyePositionTrace(t-1) - saccadeTargets(numCompletedSaccades+1);
            timeStepSaccadeMagnitude = (dt - reduceSaccadeTimeInPresentTimeStepWith)*saccadeSpeed;
            
            if (timeStepSaccadeMagnitude > saccadeOffset),
                
                % yes, so lets not do the whole thing
                eyePositionTrace(t) = saccadeTargets(numCompletedSaccades+1);
                
                % and lets saccade as completed
                numCompletedSaccades = numCompletedSaccades + 1;
                
            else
                
                % no, then lets do the full thing and keep going.
                eyePositionTrace(t) = eyePositionTrace(t-1) + -1*sign(saccadeOffset)*saccadeSpeed*dt;
                
            end
            
        end
        
    end

    % Generate retinal target location trace
    retinalTargetTraces = zeros(maxNumberOfVisibleTargets, numTimeSteps);
    for h=1:maxNumberOfVisibleTargets,
        
        % Make trace: r = h - e
        retinalTargetTraces(h, :) = headCenteredTargetLocations(h) - eyePositionTrace;
        
        % Cancel out parts where target is not present
        offIntervals = targetOffIntervals{h};
        [numIntervals x] = size(offIntervals);
        for i=1:numIntervals,
            
            % Get interval
            interval = offIntervals(i,:);
            
            % Translate from time to timesteps
            timeStepInterval = timeToTimeStep(interval);
            
            % Cancel out, i.e. not visible
            retinalTargetTraces(h, timeStepInterval(1):timeStepInterval(2)) = nan;
        end
    end
    
    function r = stepToTime(i)
        r = (i-1)*dt;
    end
    
    function i = timeToTimeStep(t)
        i = floor(t/dt) + 1;
    end

    figure;
    plot(eyePositionTrace, 'r');
    hold on;
    plot(retinalTargetTraces' , 'b');
    xlabel('Time step');
    legend({'Eye Position','Stimuli Retinal Locations'})
    %ylim(min(min(eyePositionTrace),min(min(retinalTargetTraces)) max(max(eyePositionTrace),max(max(retinalTargetTraces)))]);
    
    % Save params
    stimuliFolder = [base 'Stimuli' filesep Name];
    mkdir(stimuliFolder);
    save([stimuliFolder filesep 'TrainingStimuli.mat'] , ...
                                    'headCenteredTargetLocations', ...
                                    'maxNumberOfVisibleTargets', ...
                                    'targetOffIntervals', ...
                                    'initialEyePosition', ...
                                    'saccadeSpeed', ...
                                    'saccadeTimes', ...
                                    'saccadeTargets', ...
                                    'numSaccades', ...
                                    'eyePositionTrace', ...
                                    'retinalTargetTraces', ...
                                    'Duration', ...
                                    'dt', ...
                                    'numTimeSteps', ...
                                    'seed');
end