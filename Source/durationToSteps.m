
%
%  durationToSteps.m
%  Remapping
%
%  Created by Bedeho Mender on 11/09/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function i = durationToSteps(t,dt)
    i = floor(t/dt);
end