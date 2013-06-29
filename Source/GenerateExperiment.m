
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
    global STIMULI_FOLDER;
    
    % Stimuli
    trainingStimuli     = 'basic-Training';
    testing_kusonoki    = 'basic-KusonokiTesting';
    testing_stim_ctrl   = 'basic-StimuliControlTask';
    testing_sac_ctrl    = 'basic-SaccadeControlTask';
    
    trainingStimuliFile         = [STIMULI_FOLDER trainingStimuli filesep 'stim.mat'];
    testingKusonokiStimuliFile  = [STIMULI_FOLDER testing_kusonoki filesep 'stim.mat'];
    testingStimCTRLStimuliFile  = [STIMULI_FOLDER testing_stim_ctrl filesep 'stim.mat'];
    testingSacCTRLStimuliFile   = [STIMULI_FOLDER testing_sac_ctrl filesep 'stim.mat'];
    
    % Experiment parameters
    Name = 'prewired';
    experimentFolderPath = [EXPERIMENTS_FOLDER Name];
    
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
    dt = 0.010; % (s)
    numTrainingEpochs = 1;
    outputSavingRate = 1; % Period of time step saving during testing.
    saveActivityInTraining = false;
    saveNetworksAtEpochMultiples = 333; % Save network at this resolution
    seed = 13;
    
    % R
    parameterCombinations('R_eccentricity') = [45];
    parameterCombinations('R_tau')          = [0.100]; % (s)
    parameterCombinations('R_w_INHB')       = [0]; %0.7
    parameterCombinations('R_slope')        = [1];
    parameterCombinations('R_threshold')    = [2.0];
    parameterCombinations('R_to_C_alpha')   = [0.1]; % learning rate
    
    % V
    parameterCombinations('V_sigma')        = [5]; % (deg) receptive field size
    parameterCombinations('V_tau')          = [0.400]; % (s)
    parameterCombinations('V_psi')          = [4];
    parameterCombinations('V_to_R_psi')     = [4];
    parameterCombinations('V_to_C_psi')     = [1.0]; % R_to_C_psi
    
    % S
    parameterCombinations('S_eccentricity') = [30];
    parameterCombinations('S_delay_sigma')  = [0.4]; % (s)
    parameterCombinations('S_tau')          = [0.300]; % (s)
    parameterCombinations('S_psi')          = [1];
    parameterCombinations('S_sigma')        = parameterCombinations('V_sigma'); % (deg) receptive field size
    parameterCombinations('S_slope')        = [6];
    parameterCombinations('S_threshold')    = [0.2];
    parameterCombinations('S_to_C_psi')     = [6.2];
    parameterCombinations('S_to_C_alpha')   = [0.1]; % learning rate
    
    % C
    parameterCombinations('C_tau')          = [0.100]; % (s)
    parameterCombinations('C_w_INHB')       = [0]; %/C_N
    parameterCombinations('C_slope')        = [50];
    parameterCombinations('C_threshold')    = [1.0];
    parameterCombinations('C_to_R_psi')     = [0.15]; % 0.15
    parameterCombinations('C_to_R_alpha')   = [0.1]; % learning rate
    
    % Save the experiment params
    save([experimentFolderPath filesep 'GenerateExperiment.mat'], 'parameterCombinations');
    
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
            
            disp(['Making simulation: ' simulationName]);
            
            % Create simulation folder
            simulationFolder = [experimentFolderPath filesep simulationName];
            mkdir(simulationFolder);
            
            % Derive new paramters
            simulation('R_preferences')         = -simulation('R_eccentricity'):1:simulation('R_eccentricity');
            simulation('S_preferences')         = -simulation('S_eccentricity'):1:simulation('S_eccentricity');
            offset                              = simulation('S_delay_sigma')*randn(1, length(simulation('S_preferences'))); % (s)
            offset(offset < 0)                  = -offset(offset < 0);
            simulation('S_presaccadicOffset')   = offset;

            % Save parameters, add miscelanous paramters
            parameterfile = [simulationFolder filesep 'Parameters.mat'];
            save(parameterfile, 'simulation', 'dt', 'numTrainingEpochs', 'outputSavingRate', 'saveActivityInTraining', 'saveNetworksAtEpochMultiples', 'seed');
            
            %{
            % Create simulation blank network
            disp('Creating blank network...');
            createBlankNetwork([simulationFolder filesep 'BlankNetwork.mat'], simulation);
                        
            % Training
            disp('Training...');
            Remapping(simulationFolder, trainingStimuliFile, true);
            %}
            
            % Create prewired network
            disp('Create prewired network...');
            CreatePrewiredNetwork([simulationFolder filesep 'PrewiredNetwork.mat'], simulation('R_preferences'), simulation('S_preferences'), simulation('V_sigma'));
            
            % Move each network to new folder & test
            listing = dir(simulationFolder); 
            for d = 1:length(listing),
                
                % We looking for networks files
                subsim_name = listing(d).name;
                
                if ~listing(d).isdir && ~isempty(findstr(subsim_name,'Network')),
                    
                    % Make dir name and dir
                    networkfile = [simulationFolder filesep subsim_name];
                    [pathstr, name, ext] = fileparts(networkfile);
                    subsim_dir = [simulationFolder filesep name];
                    mkdir(subsim_dir);
                    
                    % Move files into dir
                    movefile(networkfile, subsim_dir);
                    copyfile(parameterfile, subsim_dir)
                    
                    % Testing network
                    disp('Kusonoki Task...');
                    Remapping(subsim_dir, testingKusonokiStimuliFile, false, 'kusonoki', [name ext]);
                    
                    disp('Saccade Control Task...');
                    Remapping(subsim_dir, testingSacCTRLStimuliFile, false, 'saccade-control', [name ext]);
                    
                    disp('Stimulus Control Task...');
                    Remapping(subsim_dir, testingStimCTRLStimuliFile, false, 'stimulus-control', [name ext]);

                end
            end
        end
    end
    
end