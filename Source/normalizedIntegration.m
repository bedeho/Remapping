
%
%  normalizedIntegration.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function response = normalizedIntegration(activity, dt, startTime, Duration)

    timeSteps           = timeToTimeStep(startTime, dt):1:timeToTimeStep(startTime+Duration, dt);
    offset_activity     = activity(timeSteps);
    offset_response     = squeeze(trapz(offset_activity,2)); % Integrate to find response
    response            = offset_response/(length(timeSteps) - 1); % Normaliztion step, gives normalized (sp/s) units to response
    
end