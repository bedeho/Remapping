
%
%  Sprattling_CreateBlankNetwork.m
%  Remapping
%
%  Created by Bedeho Mender on 20/08/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function C_N = Sprattling_CreateBlankNetwork(outputfile, R_preferences, C_N, R_N, S_N, R_to_R_pos_sigma, R_to_R_neg_sigma, S_to_C_connectivity, V_to_C_connectivity, C_to_R_connectivity, K_to_R_connectivity)

    %C_N = ceil(S_N*R_N*0.1); % *0.1
    
    % easy starting point
    K_N = R_N;

    % Dilution map
    C_to_R_weights_dilutionmap = createDilutionMap(R_N, C_N, C_to_R_connectivity);
    S_to_C_weights_dilutionmap = createDilutionMap(C_N, S_N, S_to_C_connectivity);
    V_to_C_weights_dilutionmap = createDilutionMap(C_N, R_N, V_to_C_connectivity);
    
    K_to_R_weights_dilutionmap = createDilutionMap(R_N, K_N, K_to_R_connectivity);
    
    % Random
    S_to_C_weights = rand(C_N,S_N).*S_to_C_weights_dilutionmap;
    V_to_C_weights = rand(C_N,R_N).*V_to_C_weights_dilutionmap;
    K_to_R_weights = rand(R_N,K_N).*K_to_R_weights_dilutionmap;
    C_to_R_weights = rand(R_N,C_N).*C_to_R_weights_dilutionmap;    
    
    %% Normalize
    S_to_C_weights = normalizeWeightVector(S_to_C_weights);
    V_to_C_weights = normalizeWeightVector(V_to_C_weights);
    
    % separate
    %C_to_R_weights = normalizeWeightVector(C_to_R_weights);
    %K_to_R_weights = normalizeWeightVector(K_to_R_weights);
    
    % combined
    [C_to_R_weights, K_to_R_weights] = Sprattling_normalizeWeightVector(C_to_R_weights, K_to_R_weights);
    
    %% Other
    % R_to_R
    [X_2 Y_2]      = meshgrid(R_preferences, R_preferences);
    
    % R_to_R_weights
    R_to_R_excitatory_weights = exp(-((X_2 - Y_2).^2)./(2*R_to_R_pos_sigma^2));

    % MAKE INTO MEXICAN HAT
    R_to_R_inhibitory_weights = exp(-((X_2 - Y_2).^2)./(2*R_to_R_neg_sigma^2)) - 0.999;%1
    R_to_R_inhibitory_weights(R_to_R_inhibitory_weights > 0) = 0;
    
    % Save params
    save(outputfile, 'C_to_R_weights', ...
                     'S_to_C_weights', ...
                     'V_to_C_weights', ...
                     'K_to_R_weights', ...
                     'R_to_R_excitatory_weights', ...
                     'R_to_R_inhibitory_weights', ...
                     'R_N', 'S_N', 'C_N', 'K_N', ...
                     'C_to_R_weights_dilutionmap', ...
                     'S_to_C_weights_dilutionmap', ...
                     'V_to_C_weights_dilutionmap', ...
                     'K_to_R_weights_dilutionmap');
    
end