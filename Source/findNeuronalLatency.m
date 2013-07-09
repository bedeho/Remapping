
%
%  findNeuronalLatency.m
%  Remapping
%
%  Created by Bedeho Mender on 09/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [latencyTimeStep, durationTimeSteps] = findNeuronalLatency(responseThreshold, responseVector, latencyWindowLength)

    % Working vars
    responseLength      = length(responseVector);
    foundOnset          = false; % Inidicator for if we have found onset
    foundOffset         = false; % Inidicator for if we have found onset/offset
    thresholdResponse   = responseVector*responseThreshold; % cutoff

    for t=1:(responseLength-latencyWindowLength+1),

        % Get time steps in question
        latencyWindow = responseVector(t + 1:latencyWindowLength);

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
                durationTimeSteps = t;
                foundOffset = true;

                return;
            end
        end
    end
end