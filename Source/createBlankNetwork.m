
%
%  createBlankNetwork.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function createBlankNetwork(outputfile, simulation)

    R_N = length(simulation('R_preferences'));
    S_N = length(simulation('S_preferences'));
    C_N = S_N*R_N;

    % Random
    C_to_R_weights = rand(R_N,C_N);
    V_to_C_weights = rand(C_N,R_N);
    S_to_C_weights = rand(C_N,S_N);
    
    % Normalize
    C_to_R_norm = 1./sqrt(squeeze(sum(C_to_R_weights.^2))); 
    C_to_R_weights = bsxfun(@times,C_to_R_weights,C_to_R_norm);

    S_to_C_norm = 1./sqrt(squeeze(sum(S_to_C_weights.^2))); 
    S_to_C_weights = bsxfun(@times,S_to_C_weights,S_to_C_norm); 

    V_to_C_norm = 1./sqrt(squeeze(sum(V_to_C_weights.^2))); 
    V_to_C_weights = bsxfun(@times,V_to_C_weights,V_to_C_norm);
    
    % Save params
    save(outputfile , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights');
    
end