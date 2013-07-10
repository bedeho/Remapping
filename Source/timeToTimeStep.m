
%
%  timeToTimeStep.m
%  Remapping
%
%  Created by Bedeho Mender on 30/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function i = timeToTimeStep(t,dt)
    i = floor(t/dt) + 1;
end