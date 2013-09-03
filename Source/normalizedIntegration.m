
%
%  normalizedIntegration.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function response = normalizedIntegration(activity, dt, startTime, Duration)

    timeSteps                 = timeToTimeStep(startTime, dt):1:timeToTimeStep(startTime+Duration, dt);
    [numNeurons numTimeSteps] = size(activity);

    if(numNeurons > 1),
        offset_activity = activity(:,timeSteps);
    else
        offset_activity = activity(timeSteps);
    end
    
    offset_response     = dt*squeeze(trapz(offset_activity,2)); % Integrate to find response
    duration            = dt*(length(timeSteps) - 1);
    response            = offset_response/duration; % Normaliztion step, gives normalized (sp/s) units to response
    
end