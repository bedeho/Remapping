
%
%  normalizeWeightVector.m
%  Remapping
%
%  Created by Bedeho Mender on 26/10/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function w = normalizeWeightVector(weightVector)

    %% ROWS MUST BE POSTSYNAPTIC NEURONS
    [num_postsynaptic num_presynaptic] = size(weightVector);

    norm = 1./sqrt(squeeze(sum(weightVector.^2, 2)));
    w = weightVector.*repmat(norm,1,num_presynaptic);
    
end