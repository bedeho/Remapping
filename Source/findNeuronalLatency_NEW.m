
%
%  findNeuronalLatency_NEW.m
%  Remapping
%
%  Created by Bedeho Mender on 10/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function latencyTimeStep = findNeuronalLatency_NEW(responseVector, dt)

    % analysis variables
    windowDuration = 0.030;
    minResponse = 0.1; %classic: 0.1
    minDiff = 0.002; % classic: 0.01 <------------ dependson dt !!!!!

    % working variables
    windowLength = timeToTimeStep(windowDuration, dt); % time over which to demand positive slop
    differential = diff(responseVector);
    responseLength = length(differential);
    
    % default value
    latencyTimeStep = nan;
    
    % check that there is some response: we could be more sophisticated
    % here, but who cares.
    if(all(responseVector < minResponse)),
        return;
    end

    % detect match: naive
    for t=1:(responseLength-windowLength),
        
        % CLASSIC
        activity = differential(t:(t+windowLength));
        
        % are all slopes in window positive ? if so we are done
        if(all(activity > minDiff))
            latencyTimeStep = t;
            return;
        end
        
        
        %{
        % NEW
        if(normalizedIntegration(responseVector, dt, (t-1)*dt, windowDuration) > 0.002)
            latencyTimeStep = t;
            return;            
        end
        %}
    end

end