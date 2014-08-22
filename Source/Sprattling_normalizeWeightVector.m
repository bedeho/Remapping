
%
%  normalizeWeightVector.m
%  Remapping
%
%  Created by Bedeho Mender on 26/10/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

%% ROWS MUST BE POSTSYNAPTIC NEURONS
% Assume that for both inputs, the number of post synaptic neurons, that is
% rows, is the same!!

function [w1,w2] = Sprattling_normalizeWeightVector(weightVector1, weightVector2)

    % Get dimensions
    [num_postsynaptic_1 num_presynaptic_1] = size(weightVector1);
    [num_postsynaptic_2 num_presynaptic_2] = size(weightVector2);
    
    assert(num_postsynaptic_1 == num_postsynaptic_2, 'BEDEHO: Incompatible number of post-synaptic output neurons.');
    
    % Find norm across total weight vector of post synaptic neurons
    global_norm = 1./sqrt(squeeze(sum(weightVector1.^2, 2) + sum(weightVector2.^2, 2)));
    
    % Multiply through on both partial weight vectors
    w1 = weightVector1.*repmat(global_norm,1,num_presynaptic_1);
    w2 = weightVector2.*repmat(global_norm,1,num_presynaptic_2);
    
    % Check that no silly nan stuff happened (can indeed happen if all
    % weights are zero)
    assert(nnz(isnan(w1)) == 0, 'BEDEHO: Weight vector has NaN due to normalizaiton error.');
    
end