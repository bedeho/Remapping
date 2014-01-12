
%
%  findNeuronalLatency_NEW.m
%  Remapping
%
%  Created by Bedeho Mender on 10/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function latencyTimeStep = findNeuronalLatency_NEW(responseVector, dt)

    % analysis variables
    windowDuration = 0.02;
    minResponse = 0.1; % do we even need this
    minDiff = 0.01;

    % working variables
    windowLength = timeToTimeStep(windowDuration, dt); % time over which to demand positive slop
    differential = diff(responseVector);
    responseLength = length(differential);
    
    % default value
    latencyTimeStep = nan;
    
    % check that there is some response: we could be more sophisticated
    % here, but who cares.
    if(all(responseLength < minResponse)),
        return;
    end

    % detect match: naive
    for t=1:(responseLength-windowLength),
        
        % are all slopes in window positive ? if so we are done
        if(all(differential(t:(t+windowLength)) > minDiff))
            latencyTimeStep = t;
            return;
        end
    end

end