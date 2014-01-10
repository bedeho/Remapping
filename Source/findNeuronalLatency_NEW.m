
%
%  findNeuronalLatency_NEW.m
%  Remapping
%
%  Created by Bedeho Mender on 10/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function latencyTimeStep = findNeuronalLatency_NEW(responseVector, dt)

    % working variables
    responseLength = length(responseVector)
    windowLength = timeToTimeStep(0.02, dt); % time over which to demand positive slop
    differential = diff(responseVector);

    % detect match: naive
    for t=1:(responseLength-windowLength+1),
        
        % are all slopes in window positive ? if so we are done
        if(all(differential(t:(t+windowLength)) > 0)
            latencyTimeStep = t;
            return;
        end
    end
    
    latencyTimeStep = nan;

end