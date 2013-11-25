
%
%  dtRoundUpPeriod.m
%  Remapping
%
%  Created by Bedeho Mender on 23/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function i = dtRoundUpPeriod(duration,dt)
    i = ceil(duration/dt)*dt;
end