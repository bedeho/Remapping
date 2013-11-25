
%
%  GenerateExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function GenerateExperiment(Name, dt, stimulinames, trainingStimuli)

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    % Experiment parameters
    %Name = 'prewired'; 
    
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
    numTrainingEpochs = 1;%10;
    doTrain = (nargin == 4) && (numTrainingEpochs > 0);
    outputSavingRate = 1; % Period of time step saving during testing.
    assert(outputSavingRate == 1, 'outputSavingRate is not 1, all further analysis will fail');
    
    saveActivityInTraining = true;
    saveNetworksAtEpochMultiples = 333; % Save network at this resolution
    seed = 13;
    
    rng(seed);
    
    % R
    parameterCombinations('R_eccentricity')         = [45];
    parameterCombinations('R_tau')                  = [0.020]; % NEW=0.050, (s), 0.050 0.100
    parameterCombinations('R_w_INHB')               = [0.6]; %WORKS: 1.8 , classic = 0.6 NEW=0.1, 20/91 15/91 10/91 ,prewwired=0 works = 5/91,20/91
    parameterCombinations('R_slope')                = [0.5]; %CHANGE BACK: 0.002 prewired=2, classic = 0.4, 0.0005
    parameterCombinations('R_threshold')            = [3]; %CHANGE BACK: 700, prewired 0.6
    parameterCombinations('R_psi')                  = [0]; % CHANGE BACK: 2000, prewired = 0.5
    %parameterCombinations('R_covariance_threshold') = [0]; % 0.2 0.3 0.4 0.8
    
    % K
    parameterCombinations('K_tau')                  = [0.020];% 0.500
    parameterCombinations('K_onset_delay_sigma')    = [0.05];
    parameterCombinations('K_I_psi')                = [8]; %NEW:4, CHANGE BACK: 1000, 16
    parameterCombinations('K_supression_delay')     = [0]; %.05,D= Duhamle suggests even shorter, an effect begins quite early if you look PSTH
    parameterCombinations('P_tau')                  = [0.300];
    parameterCombinations('P_psi')                  = [0.8];
    K_max_onset_delay = 0.080; %(s), 0.010
    
    % V
    DELAY = 0.300
    
    parameterCombinations('V_sigma')                = [3]; % (deg) receptive field size
    parameterCombinations('V_supression_delay')     = [DELAY]; % NEW:0.200
    parameterCombinations('V_to_C_psi')             = [10]; % prewired=1
    parameterCombinations('V_to_C_alpha')           = [0.1]; %WORKS: 0.5 ____0.5, NEW=0.001
    parameterCombinations('V_to_C_connectivity')    = [0.20]; % WORKS

    % S
    parameterCombinations('S_eccentricity')         = [30];
    parameterCombinations('S_delay_sigma')          = [0.100]; % (s)
    parameterCombinations('S_tau')                  = [0.020]; % (s)
    parameterCombinations('S_psi')                  = [1];
    parameterCombinations('S_presaccadic_onset')    = [0.050]; %0.050,   0.100 is classic
    parameterCombinations('S_trace_length')         = [DELAY-0.020]; % stop saccade sooner so that you dont get FRF imprinted in V->C weights due to S and C delay and V speed
    parameterCombinations('S_to_C_psi')             = [8];
    parameterCombinations('S_to_C_alpha')           = [0.1]; %WORKS: 0.5 _____0.2 NEW=0.001 %0.003, learning rate
    parameterCombinations('S_to_C_connectivity')    = [0.20]; % WORKS 0.05
    
    % C
    parameterCombinations('C_N')                    = [1000]; % 5000 full
    parameterCombinations('C_tau')                  = [0.020]; % (s)
    parameterCombinations('C_w_INHB')               = [0.1]; %1.0 WORKS , 2.8 worked, NEW: 10/5000, classic=[10/5000, 10/5000 50/5000 100/5000  C_N = 5400
    parameterCombinations('C_threshold')            = [15]; % 15 WORKS.
    parameterCombinations('C_threshold_sigma')      = [0]; % dont change!
    parameterCombinations('C_slope')                = [100]; % prewired 100, classic= 500
    parameterCombinations('C_to_R_psi')             = [8]; %10 8 6 4, NEW: 0.09 or? 0.19, CHANGE BACK: 30, 15, prewwired=0.05,0.1 works well, classic: 0.4
    parameterCombinations('C_to_R_alpha')           = [0.1]; % WORKS (no covariance): 0.20 0.1 to small? learning rate
    parameterCombinations('C_to_R_connectivity')    = [1]; %0.5 [0 1]
    
    % Save the experiment params
    save([experimentFolderPath filesep 'GenerateExperiment.mat'], 'parameterCombinations');
    
    % Start paramters permutation
    simulation = containers.Map;
    nameComponents = cell(1, length(parameterCombinations.keys));
    valueComponents = cell(1, length(parameterCombinations.keys));
    permute(1);

    % Safety check
    dt_coeff = 10;
    assert(min(parameterCombinations('R_tau')) >= dt_coeff*dt && ...
           min(parameterCombinations('S_tau')) >= dt_coeff*dt && ...
           min(parameterCombinations('C_tau')) >= dt_coeff*dt, 'One of the time constants are to small.');
    
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
            %S_delay = randn(1, S_N);
            %S_delay = 0.06 + simulation('S_delay_sigma')*S_delay; % change mean and std
            %S_delay(S_delay < 0) = -S_delay(S_delay < 0); % flip negative delays to be positive
            S_delay = simulation('S_delay_sigma')*ones(1, S_N);
            simulation('S_presaccadicOffset') = S_delay;
            
            % K delays
            K_onset_delays = randn(1, R_N); % Sample normal distribution
            K_onset_delays(K_onset_delays < 0) = -K_onset_delays(K_onset_delays < 0); % flip negative delays to be positive
            K_onset_delays = 0.00 + simulation('K_onset_delay_sigma')*K_onset_delays; % change mean and std
            K_onset_delays(K_onset_delays > K_max_onset_delay) = K_max_onset_delay; % clip delays that are to long
            simulation('K_onset_delays') = K_onset_delays; 
            
            %% Save parameters, add miscelanous paramters
            parameterfile = [simulationFolder filesep 'Parameters.mat'];
            save(parameterfile, 'simulation', 'dt', 'numTrainingEpochs', 'outputSavingRate', 'saveActivityInTraining', 'saveNetworksAtEpochMultiples', 'seed');
            
            % Create prewired network
            disp('Create prewired network...');
            
            hardwired_pref_R = simulation('R_preferences'); % 0*ones(1,R_N);, simulation('R_preferences')
            hardwired_pref_S = simulation('S_preferences'); % 18*ones(1,S_N);, simulation('S_preferences')
            
            C_N = simulation('C_N');
            
            C_to_R_sigma = simulation('V_sigma');
            V_to_C_sigma = simulation('V_sigma');
            S_to_C_sigma = simulation('V_sigma');
            R_to_R_pos_sigma = simulation('V_sigma');
            R_to_R_neg_sigma = 10*simulation('V_sigma');
            
            % Create prewired network
            %disp('Creating prewired network...');
            %CreatePrewiredNetwork([simulationFolder filesep 'PrewiredNetwork.mat'], hardwired_pref_R, hardwired_pref_S, C_N, C_to_R_sigma, V_to_C_sigma, S_to_C_sigma, R_to_R_sigma, simulation('S_to_C_connectivity'), simulation('V_to_C_connectivity'), simulation('C_to_R_connectivity'));
            
            % Create simulation blank network
            disp('Creating blank network...');
            CreateBlankNetwork([simulationFolder filesep 'BlankNetwork.mat'], hardwired_pref_R, C_N, length(simulation('R_preferences')), length(simulation('S_preferences')), R_to_R_pos_sigma, R_to_R_neg_sigma, simulation('S_to_C_connectivity'), simulation('V_to_C_connectivity'), simulation('C_to_R_connectivity'));
            
            if doTrain,
            
                % Training
                disp('Training...');
                Remapping(simulationFolder, trainingStimuli, true);            
            end
            
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
                    for s=1:length(stimulinames),
                        
                        disp(['Testing Stimuli: ' stimulinames{s}]);
                        Remapping(subsim_dir, stimulinames{s}, false, [name ext]);                       
                    end

                end
            end
        end
    end

    AnalyzeExperiment(Name, stimulinames, trainingStimuli);
    
end