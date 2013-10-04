
%
%  createPrewiredNetwork.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function CreatePrewiredNetwork(outputfile, R_preferences, S_preferences, C_to_R_sigma, V_to_C_sigma, S_to_C_sigma, R_to_R_sigma)

    R_N = length(R_preferences);
    S_N = length(S_preferences);
    C_N = S_N*R_N;

    % Prewired
    [X Y Z] = meshgrid(R_preferences, S_preferences, R_preferences);
    
    % C_to_R_weights
    C_to_R_raw      = exp(-(((Z - (X-Y))).^2)./(2*C_to_R_sigma^2));
    C_to_R_reshaped = reshape(C_to_R_raw, C_N, R_N)';
    C_to_R_norm     = 1./sqrt(squeeze(sum(C_to_R_reshaped.^2))); 
    C_to_R_weights  = bsxfun(@times, C_to_R_reshaped, C_to_R_norm); % Normalize
    
    % V_to_C_weights
    V_to_C_raw      = exp(-((Z - X).^2)./(2*V_to_C_sigma^2));
    V_to_C_reshaped = reshape(V_to_C_raw, C_N, R_N);
    V_to_C_norm     = 1./sqrt(squeeze(sum(V_to_C_reshaped.^2))); 
    V_to_C_weights  = bsxfun(@times, V_to_C_reshaped, V_to_C_norm); % Normalize
    
    % NOT LONGER IN USE: R_to_C_weights
    R_to_C_weights  = V_to_C_weights;
    
    % ADD SOMETHING HERE TO SAVE IN A FILE WHAT PAIRS OF RETINAL AND SCCADE PAIRS EACH C NEURON CORRESPONDS TO!!
    
    % NOT LONGER IN USE: V_to_R_weights
    %{
    [X_2 Y_2]       = meshgrid(R_preferences, R_preferences);
    V_to_R_raw      = exp(-((X_2 - Y_2).^2)./(2*V_sigma^2));
    V_to_R_norm     = 1./sqrt(squeeze(sum(V_to_R_raw.^2))); 
    V_to_R_weights  = bsxfun(@times, V_to_R_raw, V_to_R_norm); % Normalize
    %}
    
    % R_to_R_weights
    [X_2 Y_2]       = meshgrid(R_preferences, R_preferences);
    gauss           = exp(-((X_2 - Y_2).^2)./(2*R_to_R_sigma^2));

    % MAKE INTO MEXICAN HAT
    R_to_R_excitatory_weights = gauss;
    R_to_R_inhibitory_weights = gauss - 1;
    R_to_R_weights  = gauss;
    %R_to_R_weights  = gauss*3 - 0.1;
    
    % S_to_C_weights
    [X Y Z]         = meshgrid(R_preferences, S_preferences, S_preferences);
    
    S_to_C_raw      = exp(-((Z - Y).^2))./(2*S_to_C_sigma^2);
    S_to_C_reshaped = reshape(S_to_C_raw, C_N, S_N);
    S_to_C_norm     = 1./sqrt(squeeze(sum(S_to_C_reshaped.^2))); 
    S_to_C_weights  = bsxfun(@times, S_to_C_reshaped, S_to_C_norm); % Normalize
    
    % Save params
    save(outputfile, 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights', 'R_to_C_weights', 'R_to_R_weights', 'R_to_R_excitatory_weights', 'R_to_R_inhibitory_weights', 'R_N', 'S_N', 'C_N'); % 'V_to_R_weights'
    
end