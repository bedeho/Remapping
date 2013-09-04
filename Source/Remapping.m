
%
%  Remapping.m
%  Remapping
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Remapping(simulationFolder, stimuliName, isTraining, networkfilename)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;

    % Parse input args
    if nargin < 4,
        networkfilename = 'BlankNetwork.mat';
    end

    % Load input files
    disp('Loading input files...');
    parameters = load([simulationFolder filesep 'Parameters.mat']);
    network  = load([simulationFolder filesep networkfilename]);
    stimuli  = load([STIMULI_FOLDER stimuliName filesep 'stim.mat']);

    % Validate input
    %assert(paramter.simulation('R_eccentricity') ~=
    %stimuli.R_eccentricity, 'R_eccentricity not identical in stimuli and
    %parameters');
    
    assert(parameters.dt == stimuli.dt, 'dt not identical in stimuli and model.');
    
    %% Simulation parameters
    dt = parameters.dt; % (s)
    rng(parameters.seed);
    
    % Num epochs
    if isTraining,
        numEpochs = parameters.numTrainingEpochs;
    else
        numEpochs = 1;
    end
    
    numPeriods = length(stimuli.stimuli);
    timeStepsInPeriod = zeros(1,numPeriods);
    maxNumberOfVisibleTargets = length(stimuli.stimuli{1}.headCenteredTargetLocations); % get from first period
    
    for period=1:numPeriods,
        timeStepsInPeriod(period) = length(stimuli.stimuli{period}.eyePositionTrace);
        
        % could make allocation conditional below, but perhaps that messes up some branch prediction stuff in matlab, better off just being stingy
        if(length(stimuli.stimuli{period}.headCenteredTargetLocations) ~= maxNumberOfVisibleTargets)
            error('Number maximally visible targets must not vary across periods, this is so model can be coded to avoid slowness by avoiding new reallocations of memory'); 
        end
    end
    
    %numTimeStepsPerEpoch = sum(timeStepsInPeriod);
    maxTimeStepsPerEpoch = max(timeStepsInPeriod); % will likely never even vary...
    outputSavingRate = parameters.outputSavingRate;
    
    % Allocate buffers for activity
    if (~isTraining || isTraining && parameters.saveActivityInTraining)
        numSavedTimeSteps = ceil(maxTimeStepsPerEpoch/outputSavingRate);
    else
        numSavedTimeSteps = 0;
    end
                                     
    %% Load network
    C_to_R_weights = network.C_to_R_weights;
    S_to_C_weights = network.S_to_C_weights;
    V_to_C_weights = network.V_to_C_weights;
    %V_to_R_weights = network.V_to_R_weights;
    R_to_R_weights = network.R_to_R_weights;
    R_to_R_excitatory_weights = network.R_to_R_excitatory_weights;
    R_to_R_inhibitory_weights = network.R_to_R_inhibitory_weights;
    

    %% Load dynamical parameters
    
    % R  =======================================
    R_eccentricity  = parameters.simulation('R_eccentricity');
    R_preferences   = parameters.simulation('R_preferences');
    R_N             = length(R_preferences);
    R_activation    = zeros(1,R_N);
    R_firingrate    = zeros(1,R_N);
    
    R_tau           = parameters.simulation('R_tau');
    R_w_INHB        = parameters.simulation('R_w_INHB');
    R_slope         = parameters.simulation('R_slope');
    R_threshold     = parameters.simulation('R_threshold');
    %R_to_C_alpha    = parameters.simulation('R_to_C_alpha'); % learning rate
    %R_to_C_psi      = parameters.simulation('R_to_C_psi'); % 4
    R_psi           = parameters.simulation('R_psi');
    R_attractor_psi = parameters.simulation('R_attractor_psi');
    R_background    = parameters.simulation('R_background');
    
    R_tau_rise      = parameters.simulation('R_tau_rise');
    R_tau_decay     = parameters.simulation('R_tau_decay');
    R_tau_sigma     = parameters.simulation('R_tau_sigma');
    R_tau_threshold = parameters.simulation('R_tau_threshold');

    R_tau_dynamic   = zeros(1,R_N);
    
    % K  =======================================
    K_tau           = parameters.simulation('K_tau');
    K_psi           = parameters.simulation('K_psi');
    K_delays        = parameters.simulation('K_delays');
    
    if(maxNumberOfVisibleTargets > 0),
        K = zeros(maxNumberOfVisibleTargets, R_N);
    else
        K = 0;
    end

    K_delaysTimeSteps = timeToTimeStep(K_delays, dt);
    
    % E  =======================================
    
    E_sigma      = parameters.simulation('E_sigma');
    E_tau_rise   = parameters.simulation('E_tau_rise');
    E_tau_decay  = parameters.simulation('E_tau_decay');
    E_tau_sigma  = parameters.simulation('E_tau_sigma');
    E_to_V_psi   = parameters.simulation('E_to_V_psi');
    E_to_R_psi   = parameters.simulation('E_to_R_psi');
    
    E            = zeros(1,R_N);
    E_tau_dynamic= zeros(1,R_N);
    
    % V =======================================
    V_tau           = parameters.simulation('V_tau');
    %V_psi           = parameters.simulation('V_psi');
    %V_sigma         = parameters.simulation('V_sigma');
    V_to_R_psi      = parameters.simulation('V_to_R_psi');
    V_to_R_alpha    = parameters.simulation('V_to_R_alpha'); % learning rate
    V_to_C_psi      = parameters.simulation('V_to_C_psi');
    V_to_C_alpha    = parameters.simulation('V_to_C_alpha'); % learning rate
    V               = zeros(1,R_N);
    V_firingrate    = zeros(1,R_N);

    V_slope         = parameters.simulation('V_slope');
    V_threshold     = parameters.simulation('V_threshold');
    
    %V_tau_decay     = parameters.simulation('V_tau_decay');
    %V_tau_sigma     = parameters.simulation('V_tau_sigma');
    
    % S  =======================================
    S_eccentricity  = parameters.simulation('S_eccentricity');
    S_preferences   = parameters.simulation('S_preferences');
    S_N             = length(S_preferences);
    
    S_activation    = zeros(1,S_N);
    S_firingrate    = zeros(1,S_N);
    
    S_delay_sigma   = parameters.simulation('S_delay_sigma');
    S_presaccadicOffset = parameters.simulation('S_presaccadicOffset'); %
    S_tau           = parameters.simulation('S_tau'); % (s)
    S_to_C_psi      = parameters.simulation('S_to_C_psi');
    %S_psi           = parameters.simulation('S_psi');
    S_slope         = parameters.simulation('S_slope');
    S_threshold     = parameters.simulation('S_threshold');
    S_to_C_alpha    = parameters.simulation('S_to_C_alpha'); % learning rate
    S_sigma         = parameters.simulation('S_sigma');; % (deg) receptive field size

    % C  =======================================
    C_N             = S_N*R_N;
    C_activation    = zeros(1,C_N);
    C_firingrate    = zeros(1,C_N);
    
    C_tau           = parameters.simulation('C_tau'); % (s)
    C_to_R_alpha    = parameters.simulation('C_to_R_alpha');
    C_to_R_psi      = parameters.simulation('C_to_R_psi');
    C_w_INHB        = parameters.simulation('C_w_INHB');
    C_slope         = parameters.simulation('C_slope');
    C_threshold     = parameters.simulation('C_threshold');
    
    % Allocate buffer space
    E_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    V_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    R_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    S_firing_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
    C_firing_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);

    E_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    V_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    R_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    S_activation_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
    C_activation_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);
    
    % Flat buffers
    C_firing_history_flat     = zeros(C_N, numPeriods, numEpochs);
    C_activation_history_flat = zeros(C_N, numPeriods, numEpochs);
    
    %% Simulate
    totalTicID = tic;
    for epoch=1:numEpochs,
        
        if isTraining,
            disp(['Starting epoch #' num2str(epoch)]);
        end
        
        for period=1:numPeriods,
            
            periodTicID = tic;
            
            % Get number of time steps in this period
            numTimeSteps = timeStepsInPeriod(period);
            
            maxNumberOfVisibleTargets = length(stimuli.stimuli{period}.headCenteredTargetLocations);
    
            % Load Stimuli for period
            retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces;
            saccadeTimes        = stimuli.stimuli{period}.saccadeTimes;
            saccadeTargets      = stimuli.stimuli{period}.saccadeTargets;
            stimOnsetTimes      = stimuli.stimuli{period}.stimOnsetTimes;
            stimOnsetTimeSteps  = timeToTimeStep(stimOnsetTimes,dt);

            numSaccades         = length(stimuli.stimuli{period}.saccadeTimes);
            numStimOnsetTimes   = length(stimOnsetTimes);
            lastStimOnsetTime   = max(stimOnsetTimes);
            
            % Setup static working variables
            R_preference_comparison_matrix = repmat(R_preferences, maxNumberOfVisibleTargets, 1);

            stimOnset_comparison_matrix = repmat(stimOnsetTimeSteps', 1, R_N); 
            
            K_delay_comparison_matrix = repmat(K_delaysTimeSteps, numStimOnsetTimes, 1);
            
            %stimulus_presence_indicator = ~isnan(retinalTargetTraces);
            
            if(numSaccades > 0),
                saccade_times = repmat(saccadeTimes,S_N,1); % Used to get F
                saccade_time_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                saccade_time_pos_offset = saccade_times + repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                saccade_time_neg_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                saccade_target_offset = exp(-((repmat(S_preferences', 1, numSaccades) - repmat(saccadeTargets, S_N, 1)).^2)./(2*S_sigma^2)); % term multiplied by F
            end
            
            % Reset network variables
            periodSaveCounter   = 2; % First dt is automatically saved
            E                   = E*0;
            V                   = V*0;
            K                   = K*0; % must be reallocated since teh number of visiible targets can in theory change, although it never does in practice.
            R_firingrate        = R_firingrate*0;
            S_firingrate        = S_firingrate*0;
            C_firingrate        = C_firingrate*0;
            V_firingrate        = V_firingrate*0;
            
            R_activation        = R_activation*0;
            S_activation        = S_activation*0;
            C_activation        = C_activation*0;
    
            % Run period
            for t=2:numTimeSteps, % we cannot compute anything for t==1, since this is time=0, which are the initial conditions, we have not t=-dt input data to use for this

                % Turn time step into real time
                %time             = stepToTime(t, dt);
                time_precedingdt = stepToTime(t-1, dt);
                
                % Time comparison for delta impulse must be done in
                % time steps, not continous time, otherwise we it will
                % almost surely miss delta(0)
                precedingTimeStep = t-1;

                %% Activation
                
                % Keep old solution variables, they are used to compute new
                % ones in Euler scheme
                E_old = E;
                V_old = V;
                K_old = K;
                
                % Prepare response of all neurons to all visible targets
                if(~isempty(retinalTargetTraces)),
                    
                    retinalTargets = retinalTargetTraces(:,t-1); % retinal locations of stimuli in teh precedng time step
                    diff = R_preference_comparison_matrix - repmat(retinalTargets,1,R_N);
                    diff(isnan(diff)) = inf; % for nan values, make inf, so that exponantiation gives no contribution
                    gauss = exp((-(diff).^2)./(2*E_sigma^2)); % gauss has dimensions: maxNumberOfVisibleTargets X R_N
                else
                    gauss = 0;
                end
                
                flat_gauss = sum(gauss,1); % collapse stimulus dimension, so that each neuron has one driving sum of exponentials, one exponential per stimulus
                
                % E & V =======================================
                
                E_tau_dynamic = E_tau_rise + exp(-(flat_gauss.^2)/(2*E_tau_sigma^2))*(E_tau_decay-E_tau_rise);
                
                E = E_old + (dt./E_tau_dynamic).*(-E_old + flat_gauss);
                
                V = V_old + (dt/V_tau)*(-V_old + E_to_V_psi*E_old);
                
                % R =======================================
                
                % K stim onset
                if(~isempty(stimOnsetTimes)),
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    
                    K_onset_spike = K_psi*flat_gauss(yes_delta_event);
                else
                    K_onset_spike = 0;
                end
                
                % K saccade onset
                if(~isempty(saccade_times) && any(saccade_times==saccade_times)),
                    K_sacc_supression = K_psi;
                else
                    K_sacc_supression = 0;
                end
                
                % K dynamics
                K_tau_dynamic = R_tau_rise + (flat_gauss <= R_tau_threshold)*(R_tau_decay-R_tau_rise);
                K = K_old + (dt./K_tau_dynamic).*(-K + R_psi*flat_gauss) + K_onset_spike - K_sacc_supression;
                
                % R dynamics
                R_inhibition = R_w_INHB*sum(R_firingrate);
                C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate')';
                R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation - R_inhibition + K - R_background);
                
                % Old
                %{
                R_attractor = 0;%R_attractor_psi*(R_to_R_weights*R_firingrate')';
                R_inhibition = R_w_INHB*sum(R_firingrate);
                C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate')';
                
                K = -R_activation + C_to_R_excitation - R_inhibition + R_psi*flat_gauss + R_attractor - R_background;
                
                if(~isempty(stimOnsetTimes)),
                    
                    % Time comparison for delta impulse must be done in
                    % time steps, not continous time, otherwise we it will
                    % almost surely miss delta(0)
                    precedingTimeStep = t-1;
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    no_delta_event = ~yes_delta_event;

                    % Update neurons with delta event: dirac delta case
                    R_activation(yes_delta_event) = R_activation(yes_delta_event) + K_psi*flat_gauss(yes_delta_event); % 

                    % Update neurons without delta event: standard Forward
                    % Euler
                    R_activation(no_delta_event) = R_activation(no_delta_event) + (dt/R_tau)*(-R_activation(no_delta_event) + C_to_R_excitation(no_delta_event) - R_inhibition + R_psi*flat_gauss(no_delta_event) + R_attractor(no_delta_event) - R_background);
                    
                else
                    
                    % Do all neuron at ones
                    R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation - R_inhibition + R_psi*flat_gauss + R_attractor - R_background);
                end
                %}
            
                % C =======================================
                
                C_inhibition = C_w_INHB*sum(C_firingrate);
                %R_to_C_excitation = R_to_C_psi*(R_to_C_weights*R_firingrate');
                V_to_C_excitation = V_to_C_psi*(V_to_C_weights*V_firingrate');
                S_to_C_excitation = S_to_C_psi*(S_to_C_weights*S_firingrate');
                C_activation = C_activation + (dt/C_tau)*(-C_activation + V_to_C_excitation' + S_to_C_excitation' - C_inhibition); % _to_C_excitation

                % S =======================================
                if(numSaccades > 0),
                    F = (saccade_time_neg_offset <= time_precedingdt) & (time_precedingdt <= saccade_time_pos_offset); % check both conditions: y-z <= x <= y
                    %F = (saccade_time_offset <= time_precedingdt) & (time_precedingdt <= saccade_times); % check both conditions: y-z <= x <= y
                    %S_driver = S_psi*sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
                    S_driver = sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
                    S_activation = S_activation + (dt/S_tau)*(-S_activation + S_driver');
                else
                    S_activation = S_activation + (dt/S_tau)*(-S_activation); 
                end
                
                %% Weight Update
                if isTraining,

                    C_to_R_weights = C_to_R_weights + dt*C_to_R_alpha*(R_firingrate'*C_firingrate);
                    S_to_C_weights = S_to_C_weights + dt*S_to_C_alpha*(C_firingrate'*S_firingrate);
                    %R_to_C_weights = R_to_C_weights + dt*R_to_C_alpha*(C_firingrate'*R_firingrate);
                    V_to_C_weights = V_to_C_weights + dt*V_to_C_alpha*(C_firingrate'*V_firingrate);
                    %V_to_R_weights = V_to_R_weights + dt*V_to_R_alpha*(R_firingrate'*V_firingrate);

                    % Normalize
                    C_to_R_norm = 1./sqrt(squeeze(sum(C_to_R_weights.^2))); 
                    C_to_R_weights = bsxfun(@times,C_to_R_weights, C_to_R_norm);

                    S_to_C_norm = 1./sqrt(squeeze(sum(S_to_C_weights.^2))); 
                    S_to_C_weights = bsxfun(@times,S_to_C_weights, S_to_C_norm); 

                    %R_to_C_norm = 1./sqrt(squeeze(sum(R_to_C_weights.^2))); 
                    %R_to_C_weights = bsxfun(@times,R_to_C_weights,R_to_C_norm);
                    
                    V_to_C_norm = 1./sqrt(squeeze(sum(V_to_C_weights.^2))); 
                    V_to_C_weights = bsxfun(@times,V_to_C_weights, V_to_C_norm);
                    
                    %V_to_R_norm = 1./sqrt(squeeze(sum(V_to_R_weights.^2))); 
                    %V_to_R_weights = bsxfun(@times,V_to_R_weights, V_to_R_norm);
                    
                end

                %% Compute firing rates
                R_firingrate = 1./(1 + exp(-2*R_slope*(R_activation - R_threshold)));
                %S_firingrate = 1./(1 + exp(-2*S_slope*(S_activation - S_threshold)));
                S_firingrate = S_activation;
                C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_threshold)));
                V_firingrate = 1./(1 + exp(-2*C_slope*(V - V_threshold))); 

                %% Save activity                
                if (~isTraining || isTraining && parameters.saveActivityInTraining) % && mod(t, outputSavingRate) == 0,
                    
                    E_firing_history(:, periodSaveCounter, period, epoch) = K;%E;
                    V_firing_history(:, periodSaveCounter, period, epoch) = V_firingrate;
                    R_firing_history(:, periodSaveCounter, period, epoch) = R_firingrate;
                    S_firing_history(:, periodSaveCounter, period, epoch) = S_firingrate;
                    C_firing_history(:, periodSaveCounter, period, epoch) = C_firingrate;

                    E_activation_history(:, periodSaveCounter, period, epoch) = E;
                    V_activation_history(:, periodSaveCounter, period, epoch) = V;
                    R_activation_history(:, periodSaveCounter, period, epoch) = R_activation;
                    S_activation_history(:, periodSaveCounter, period, epoch) = S_activation;
                    C_activation_history(:, periodSaveCounter, period, epoch) = C_activation;
                    
                    % Count one more dt
                    periodSaveCounter = periodSaveCounter + 1;
                end 
            end
            
            % HACK to save C for C Probe Task, we dont have space to do
            % this properly
            if numSaccades > 0,

                %hacked solution to get rate during C probing task,
                %despte all data not being savable.
                firstSaccadeTime = saccadeTimes(1);
                window = 0.050;

                C_firing_history_flat(:, period, epoch) = normalizedIntegration(C_firing_history(:, :, period, epoch), dt, firstSaccadeTime - window, window);
                C_activation_history_flat(:, period, epoch) = normalizedIntegration(C_activation_history(:, :, period, epoch), dt, firstSaccadeTime - window, window);
            end
            
            periodFinishTime = toc(periodTicID);
            disp(['Finished period #' num2str(period) ' of ' num2str(numPeriods) ' in ' num2str(periodFinishTime) 's.']);
            
        end
        
        % Output intermediate network
        if mod(epoch, parameters.saveNetworksAtEpochMultiples) == 0,
            
            disp('Saving trained network to disk...');
            save([simulationFolder filesep 'TrainedNetwork_e' num2str(epoch) '.mat'] , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights', 'V_to_R_weights');
        end
    end
    
    totalFinishTime = toc(totalTicID);
    disp(['Completed ' num2str(numEpochs) ' epochs in ' num2str(totalFinishTime) 's.']);
    
    % Output final network
    if isTraining,
        disp('Saving trained network to disk...');
        save([simulationFolder filesep 'TrainedNetwork.mat'] , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights', 'V_to_R_weights', 'R_N', 'S_N', 'C_N');
    end
    
    
    if numel(C_firing_history) > 2000000/2, % dont save if bigger than 100MB
        disp('Could not save C_firing_history, TO BIG !!!!!!!!!!!!!');
        C_firing_history = [];
        C_activation_history = [];
    end
    
    
    % Output activity
    disp('Saving activity...');
    
    % SaveActivity
    save([simulationFolder filesep 'activity-' stimuliName '.mat']  , 'R_N' ...
                                                                            , 'S_N' ...
                                                                            , 'C_N' ...
                                                                            , 'numEpochs' ...
                                                                            , 'numPeriods' ...
                                                                            , 'numSavedTimeSteps' ...
                                                                            , 'outputSavingRate' ...
                                                                            , 'dt' ...
                                                                            , 'E_firing_history' ...
                                                                            , 'V_firing_history' ...
                                                                            , 'R_firing_history' ...
                                                                            , 'S_firing_history' ...
                                                                            , 'C_firing_history' ...
                                                                            , 'E_activation_history' ...
                                                                            , 'V_activation_history' ...
                                                                            , 'R_activation_history' ...
                                                                            , 'S_activation_history' ...
                                                                            , 'C_activation_history' ...
                                                                            , 'C_firing_history_flat' ...
                                                                            , 'C_activation_history_flat' ...
                                                                            , 'E');
                                                                        
    disp('Done...');

end