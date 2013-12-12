
%
%  findNeuronalLatency_OLD.m
%  Remapping
%
%  Created by Bedeho Mender on 09/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [latencyTimeStep, duration] = findNeuronalLatency(responseVector, latencyWindowLength)

    % analysis params, unifrom across all purposes
    responseThreshold = 0.4;
    minimumThreshold = 0.1;

    % Working vars
    responseLength      = length(responseVector);
    foundOnset          = false; % Inidicator for if we have found onset
    foundOffset         = false; % Inidicator for if we have found onset/offset
    thresholdResponse   = max(responseVector)*responseThreshold; % cutoff
    
    latencyTimeStep = nan;
    duration = nan;
    
    % Check for minimal discharge
    if(all(responseVector < minimumThreshold))
        return;
    end
    
    for t=1:(responseLength-latencyWindowLength+1),

        % Get time steps in question
        latencyWindow = responseVector(t + (0:(latencyWindowLength-1)));

        % Integrate to find response
        windowResponse = trapz(latencyWindow);

        % Normalize
        windowResponse = windowResponse/(latencyWindowLength-1);

        if ~foundOnset,

            if windowResponse >= thresholdResponse;
                latencyTimeStep = t;
                foundOnset = true;
            end

        elseif ~foundOffset,

            if windowResponse <= thresholdResponse;
                duration = t-latencyTimeStep;
                foundOffset = true;

                return;
            end
        end
    end

end