
%
%  GenerateTrace.m
%  Remapping
%
%  Created by Bedeho Mender on 23/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [eyePositionTrace, retinalTargetTraces] = GenerateTrace(Duration, dt, headCenteredTargetLocations, targetOffIntervals, initialEyePosition, saccadeTimes, saccadeTargets)

    numTimeSteps = timeToTimeStep(Duration,dt);
    
    % Visual
    maxNumberOfVisibleTargets = length(headCenteredTargetLocations);
    assert(length(targetOffIntervals) >= maxNumberOfVisibleTargets, 'On off history not provided for all targets.');
    
    % Saccadic
    saccadeSpeed = 300; % (s), %if changed, then change in GenerateKusonokiTestingStimuli.m as well!
    numSaccades = length(saccadeTimes);
    
    if numSaccades > 0,
        
        assert(length(saccadeTimes) >= length(saccadeTargets), 'Number of saccade times and targets must match.');
        assert(Duration >= max(saccadeTimes), 'Saccade time after total duration found.');
    end

    % Allocate space for traces
    % Index i => Time (i-1)*dt
    eyePositionTrace = zeros(1, numTimeSteps);
    eyePositionTrace(1) = initialEyePosition;

    % Generate eye position trace
    numCompletedSaccades = 0;
    for t=2:numTimeSteps,
        
        presentTime = stepToTime(t,dt);
        
        % Have we completed all saccades beginning before time timestep t?
        if numCompletedSaccades == numSaccades || presentTime <= saccadeTimes(numCompletedSaccades+1),

            % yes, so lets just continue fixating
            eyePositionTrace(t) = eyePositionTrace(t-1);

        else
            % no, there was a not completed saccade beginning before
            % timestep t
            
            % was this the first step across this saccade onset time?
            previousTimeSteptime = stepToTime(t-1,dt);
            timeOffset = abs(previousTimeSteptime - saccadeTimes(numCompletedSaccades+1));
            
            if(timeOffset < dt),
                
                % yes it was, so mini saccade starting
                % point needs to correct for delay between prior time step
                % and saccade onset event to find saccade start position
                
                reduceSaccadeTimeInPresentTimeStepWith = timeOffset;
                
                % Keep track of where we want to go
                newEyePosition = eyePositionTrace(t-1) + saccadeTargets(numCompletedSaccades+1);
            else
                
                % no, we passed at some previous time, so mini saccade starting
                % point is just last eye position
                
                reduceSaccadeTimeInPresentTimeStepWith = 0;
            end
            
            % Will one more stime step saccading from
            % miniSaccadeStartPosition take us past saccade target?
            
            remainingEyePositionOffset = newEyePosition - eyePositionTrace(t-1);
            timeStepSaccadeMagnitude = (dt - reduceSaccadeTimeInPresentTimeStepWith)*saccadeSpeed;
            
            if (timeStepSaccadeMagnitude > remainingEyePositionOffset),
                
                % yes, so lets not do the whole thing
                eyePositionTrace(t) = newEyePosition;
                
                % and lets saccade as completed
                numCompletedSaccades = numCompletedSaccades + 1;
                
            else
                
                % no, then lets do the full thing and keep going.
                eyePositionTrace(t) = eyePositionTrace(t-1) + sign(remainingEyePositionOffset)*saccadeSpeed*dt;
                
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
            timeStepInterval = timeToTimeStep(interval,dt);
            
            % Cancel out, i.e. not visible
            retinalTargetTraces(h, timeStepInterval(1):timeStepInterval(2)) = nan;
        end
    end
end