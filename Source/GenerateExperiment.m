
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
    Name = 'prewired'; 
    
    % first-training
    % kusonokireal-prewired-tuning

    % Create experiment folders
    experimentFolderPath = [EXPERIMENTS_FOLDER Name];
    
    if exist(experimentFolderPath),
        
        if strcmp(questdlg('DELETE OLD EXPERIMENT?', 'Delete', 'NO','YES','NO'), 'NO'),
            return;
        else
            rmdir(experimentFolderPath,'s')
        end   
    end
    
    mkdir(experimentFolderPath);
    
    %% Specify main paramters
    parameterCombinations = containers.Map;
    
    % Simulations Paramters
    dt = 0.010; % (s)
    numTrainingEpochs = 1;
    outputSavingRate = 1; % Period of time step saving during testing.
    assert(outputSavingRate == 1, 'outputSavingRate is not 1, all further analysis will fail');
    
    saveActivityInTraining = false;
    saveNetworksAtEpochMultiples = 333; % Save network at this resolution
    seed = 13;
    
    % R
    parameterCombinations('R_eccentricity') = [45];
    parameterCombinations('R_tau')          = [0.050]; % (s)
    parameterCombinations('R_w_INHB')       = [5/91]; %0.7 20/91 15/91 10/91 
    parameterCombinations('R_slope')        = [0.7];
    parameterCombinations('R_threshold')    = [2.0];
    %parameterCombinations('R_to_C_alpha')  = [0.1]; % learning rate
    %parameterCombinations('R_to_C_psi')    = [1];
    parameterCombinations('R_psi')          = [1.0];
    
    parameterCombinations('R_tau_rise')     = [0.050];
    parameterCombinations('R_tau_decay')    = [0.700];
    parameterCombinations('R_tau_sigma')    = [5];
    parameterCombinations('R_tau_threshold')= [0.4];
    
    % K
    parameterCombinations('K_tau')          = [0.700];
    parameterCombinations('K_psi')          = [4 5 6];
    parameterCombinations('K_delay_sigma')  = [0.05];    
    
    % E
    parameterCombinations('E_sigma')        = [5]; % (deg) receptive field size
    parameterCombinations('E_tau_rise')     = [0.05];
    parameterCombinations('E_tau_decay')    = [0.7];
    parameterCombinations('E_to_V_psi')     = [1];
    parameterCombinations('E_to_R_psi')     = [6];
    
    % V
    %parameterCombinations('V_sigma')        = [5]; % (deg) receptive field size
    parameterCombinations('V_tau')           = [0.050]; % (s)
    %parameterCombinations('V_psi')          = [1];
    %parameterCombinations('V_slope')        = [1];
    %parameterCombinations('V_threshold')    = [0.5];
    
    parameterCombinations('V_to_R_psi')     = [6]; % 5 works
    parameterCombinations('V_to_R_alpha')   = [0.1];
    
    parameterCombinations('V_to_C_psi')     = [1];
    parameterCombinations('V_to_C_alpha')   = [0.1];

    % S
    parameterCombinations('S_eccentricity') = [30];
    parameterCombinations('S_delay_sigma')  = [0.100]; % (s)
    parameterCombinations('S_tau')          = [0.020]; % (s)
    parameterCombinations('S_sigma')        = parameterCombinations('E_sigma'); % (deg) receptive field size
    %parameterCombinations('S_psi')          = [1];
    parameterCombinations('S_slope')        = [10];
    parameterCombinations('S_threshold')    = [0.5];
    
    parameterCombinations('S_to_C_psi')     = [2.3];
    parameterCombinations('S_to_C_alpha')   = [0.1]; % learning rate
    
    % C
    parameterCombinations('C_tau')          = [0.010]; % (s)
    parameterCombinations('C_w_INHB')       = [1/5000]; %10/5000 50/5000 100/5000  C_N = 5400
    parameterCombinations('C_slope')        = [1000000]; % classic= 500
    parameterCombinations('C_threshold')    = [0.5]; % old 0.45
    parameterCombinations('C_to_R_psi')     = [0.125]; % classic: 0.5
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
                    
                    if strcmp(simulationName,''),
                        simulationName = [nameComponents{p} '=' num2str(valueComponents{p})];
                    else
                        simulationName = [simulationName '-' nameComponents{p} '=' num2str(valueComponents{p})];
                    end
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
            
            %% Derive new paramters
           
            simulation('R_preferences') = -simulation('R_eccentricity'):1:simulation('R_eccentricity');
            simulation('S_preferences') = -simulation('S_eccentricity'):1:simulation('S_eccentricity');
            
            R_N = length(simulation('R_preferences'));
            S_N = length(simulation('S_preferences'));
            
            % S delays
            S_delay = randn(1, S_N);
            S_delay(S_delay < 0) = -S_delay(S_delay < 0); % flip negative delays to be positive
            S_delay = 0.1 + simulation('S_delay_sigma')*S_delay; % change mean and std
            %S_delay = ones(1, S_N)*simulation('S_delay_sigma'); %unifromity, vs. distribution = offset;
            simulation('S_presaccadicOffset') = S_delay;
            
            % K delays
            K_delays = randn(1, R_N); % Sample normal distribution
            K_delays(K_delays < 0) = -K_delays(K_delays < 0); % flip negative delays to be positive
            K_delays = 0.00 + simulation('K_delay_sigma')*K_delays; % change mean and std
            K_delays(K_delays > 0.1) = 0.1; % clip delays that are to long
            simulation('K_delays') = K_delays; 
            
            % Sigma for decays, they are derived from thresholds
            %simulation('V_tau_sigma')    = 0.5*(simulation('V_psi')*exp(-1/2)); % when a V neuron has drive
            simulation('E_tau_sigma') = 0.5*(simulation('E_to_V_psi')*exp(-1/2)); % time constant switch sigma is set to standard deviation of E tuning curve, then
            simulation('R_tau_sigma') = 0.1*(simulation('R_psi')*exp(-1/2)); % time constant switch sigma is set to standard deviation of R tuning curve, then
            
            %% Save parameters, add miscelanous paramters
            parameterfile = [simulationFolder filesep 'Parameters.mat'];
            save(parameterfile, 'simulation', 'dt', 'numTrainingEpochs', 'outputSavingRate', 'saveActivityInTraining', 'saveNetworksAtEpochMultiples', 'seed');
            
            %{
            % Create simulation blank network
            disp('Creating blank network...');
            CreateBlankNetwork([simulationFolder filesep 'BlankNetwork.mat'], length(simulation('R_preferences')), length(simulation('S_preferences')));
                        
            % Training
            disp('Training...');
            Remapping(simulationFolder, 'basic-Training', true);
            %}
            
            % Create prewired network
            disp('Create prewired network...');
            
            hardwired_pref_R = simulation('R_preferences'); % 0*ones(1,R_N);, simulation('R_preferences')
            hardwired_pref_S = simulation('S_preferences'); % 18*ones(1,S_N);, simulation('S_preferences')
            
            C_to_R_sigma = simulation('E_sigma');
            V_to_C_sigma = simulation('E_sigma');
            S_to_C_sigma = simulation('E_sigma');
            
            CreatePrewiredNetwork([simulationFolder filesep 'PrewiredNetwork.mat'], hardwired_pref_R, hardwired_pref_S, C_to_R_sigma, V_to_C_sigma, S_to_C_sigma);
            
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
                    copyfile(parameterfile, subsim_dir);
                    
                    % Testing network
                    disp('Kusonoki ...');
                    Remapping(subsim_dir, 'basic-Kusonoki', false, [name ext]);
                    
                    disp('Duhamel Remapping ...');
                    Remapping(subsim_dir, 'basic-DuhamelRemapping', false, [name ext]);
                    
                    disp('Duhamel Remapping Trace ...');
                    Remapping(subsim_dir, 'basic-DuhamelRemappingTrace', false, [name ext]);
                    
                    disp('Duhamel Truncation ...');
                    Remapping(subsim_dir, 'basic-DuhamelTruncation', false, [name ext]);
                    
                    disp('Saccade Control ...');
                    Remapping(subsim_dir, 'basic-SaccadeControl', false, [name ext]);
                    
                    disp('Stimulus Control ...');
                    Remapping(subsim_dir, 'basic-StimuliControl', false, [name ext]);
                    
                    disp('C Layer Probe ...');
                    Remapping(subsim_dir, 'basic-CLayerProbe', false, [name ext]);

                end
            end
        end
    end

    AnalyzeExperiment(Name)
    
end