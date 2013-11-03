
%
%  createDilutionMap.m
%  Remapping
%
%  Created by Bedeho Mender on 02/10/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function map = createDilutionMap(numPostSynaptic, numPresynaptic, connectivityRate)

    numAfferents = ceil(numPresynaptic*connectivityRate);

    % check that neurons will have at least one afferent, otherwise
    % weight normalization will create NaN hell
    if numAfferents < 1,
        error('Connectivity rate is to low.'); % annoying bug
    end

    % Setup map
    map = zeros(numPostSynaptic, numPresynaptic);

    % Iterate post synaptic neurons and enable connections at desired rate
    for p=1:numPostSynaptic,
        on = randperm(numPresynaptic,numAfferents);
        map(p, on) = 1;
    end
end
