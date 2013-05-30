
%
%  GenerateExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateExperiment()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;

    % Experiment parameters
    Name = 'test';
    experimentFolderPath = [EXPERIMENTS_FOLDER Name];
    
    % Stimuli
    %trainingStimuli = 'dddd';
    %testingStimuli = 'dddd';
    
    % Create experiment folders
    if exist(experimentFolderPath),
        
        if strcmp(questdlg('DELETE OLD EXPERIMENT?', 'Delete', 'NO','YES','NO'), 'NO'),
            return;
        else
            system(['rm -R ' experimentFolderPath]);
        end   
    end
    
    mkdir(experimentFolderPath);
    
    %% Specify main paramters
    parameterCombinations = containers.Map;
    
    % Simulations Paramters
    dt = 0.001; % (s)
    epochs = 1;
    outputSavingRate = 2; % Period of time step saving during testing.
    saveDuringTraining = false;
    saveNetworksAtEpochs = 333; % Save network at this resolution
    
    
    % LIP
    parameterCombinations('R_eccentricity') = [45 50];
    parameterCombinations('R_tau')          = [0.100 0.2]; % (s)
    parameterCombinations('C_to_R_psi')     = [0.08]; % 0.15
    parameterCombinations('R_w_INHB')       = [0]; %0.7
    parameterCombinations('V_tau')          = [0.400]; % (s)
    parameterCombinations('V_psi')          = [4];
    parameterCombinations('V_sigma')        = [5]; % (deg) receptive field size
    parameterCombinations('C_to_R_alpha')   = [0.1]; % learning rate
    parameterCombinations('R_slope')        = [10];
    parameterCombinations('R_threshold')    = [2.0];
    
    % FEF: Saccade Plan
    parameterCombinations('S_eccentricity') = [30];
    parameterCombinations('S_delay_sigma')  = [0.4]; % (s)
    parameterCombinations('S_tau')          = [0.300]; % (s)
    parameterCombinations('S_psi')          = [1];
    parameterCombinations('S_sigma')        = parameterCombinations('V_sigma'); % (deg) receptive field size
    parameterCombinations('S_slope')        = [6];
    parameterCombinations('S_threshold')    = [0.2];
    
    % SC?: Comb
    parameterCombinations('C_tau')          = [0.100]; % (s)
    parameterCombinations('R_to_C_psi')     = [1.0];
    parameterCombinations('S_to_C_psi')     = [6.2];
    parameterCombinations('C_w_INHB')       = [0]; %/C_N
    parameterCombinations('R_to_C_alpha')   = [0.1]; % learning rate
    parameterCombinations('S_to_C_alpha')   = [0.1]; % learning rate
    parameterCombinations('C_slope')        = [50];
    parameterCombinations('C_threshold')    = [1.0];
    
    % Save the experiment params
    save([experimentFolderPath filesep 'GenerateExperiment.mat'] , 'parameterCombinations');
    
    % Start paramters permutation
    simulation = containers.Map;
    nameComponents = cell(1, length(parameterCombinations.keys));
    valueComponents = cell(1, length(parameterCombinations.keys));
    permute(1);

    %% Generate simulations
    function permute(paramnr)
        
        keys = parameterCombinations.keys;
        
        % Are we done
        if paramnr <= length(keys),
            
            key = keys{paramnr};
            values = parameterCombinations(key);
            
            % Add to name
            if length(values) > 1,
                nameComponents{paramnr} = key;
            end
            
            for v=values,
                valueComponents{paramnr} = v;
                simulation(key) = v;
                permute(paramnr+1);
            end
        else
            
            % Make name
            simulationName = '';
            for p=1:length(nameComponents),
                
                if ~isempty(nameComponents{p}),
                    simulationName = [simulationName '-' nameComponents{p} '=' num2str(valueComponents{p})];
                end
            end
            
            % If there is only one param combination, then just name it
            % blank
            if isempty(simulationName),
                simulationName = 'baseline';
            end
            
            % Create simulation folder
            simulationFolder = [experimentFolderPath filesep simulationName];
            mkdir(simulationFolder);
            
            % Derive new paramters
            simulation('R_preferences')         = -simulation('R_eccentricity'):1:simulation('R_eccentricity');
            simulation('S_preferences')         = -simulation('S_eccentricity'):1:simulation('S_eccentricity');
            offset                              = 0.4*randn(1, length(simulation('S_preferences'))); % (s)
            offset(offset < 0)                  = -offset(offset < 0);
            simulation('S_presaccadicOffset')   = offset;
            
            % Save params
            save([simulationFolder filesep 'SimulationExperiment.mat'] , 'simulation');
            
            % Create simulation blank network
            createBlankNetwork([simulationFolder filesep 'BlankNetwork.mat'], simulation);
            
            % Training
            %Remapping(sourcefolder, outputFolder, paramters,enablePlasticity)
            
            % Testing
            %Remapping(sourcefolder, outputFolder, paramters,enablePlasticity)
            
            % Analysis
            
        end
    end
    
end