
%
%  Remapping.m
%  Remapping
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Remapping(simulationFolder, stimuliFile, isTraining, outputpostfix, networkfilename)

    % Parse input args
    if nargin < 5,
        networkfilename = 'BlankNetwork.mat';

        if nargin < 4,
            outputpostfix = '';
        else
            outputpostfix = ['-' outputpostfix];
        end
    end

    % Load input files
    disp('Loading input files...');
    parameters = load([simulationFolder filesep 'Parameters.mat']);
    network  = load([simulationFolder filesep networkfilename]);
    stimuli  = load(stimuliFile);

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
    
    for period=1:numPeriods,
        timeStepsInPeriod(period) = length(stimuli.stimuli{period}.eyePositionTrace);
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

    %% Load dynamical parameters
    
    % LIP: Retina
    R_eccentricity  = parameters.simulation('R_eccentricity');
    R_preferences   = parameters.simulation('R_preferences');
    R_N             = length(R_preferences);
    R_activation    = zeros(1,R_N);
    R_firingrate    = zeros(1,R_N);
    
    R_tau           = parameters.simulation('R_tau');
    R_w_INHB        = parameters.simulation('R_w_INHB');
    R_slope         = parameters.simulation('R_slope');
    R_threshold     = parameters.simulation('R_threshold');
    R_to_C_alpha    = parameters.simulation('R_to_C_alpha'); % learning rate
    %R_to_C_psi     = parameters.simulation('R_to_C_psi'); % 4
    
    % V
    V_tau           = parameters.simulation('V_tau');
    V_psi           = parameters.simulation('V_psi');
    V_sigma         = parameters.simulation('V_sigma');
    V               = zeros(1,R_N);
    
    V_to_R_psi      = parameters.simulation('V_to_R_psi');
    V_to_C_psi      = parameters.simulation('V_to_C_psi');

    % S
    S_eccentricity  = parameters.simulation('S_eccentricity');
    S_preferences   = parameters.simulation('S_preferences');
    S_N             = length(S_preferences);
    
    S_activation    = zeros(1,S_N);
    S_firingrate    = zeros(1,S_N);
    
    S_delay_sigma   = parameters.simulation('S_delay_sigma');
    S_presaccadicOffset = parameters.simulation('S_presaccadicOffset');
    S_tau           = parameters.simulation('S_tau'); % (s)
    S_psi           = parameters.simulation('S_psi');
    S_to_C_psi      = parameters.simulation('S_to_C_psi');
    S_slope         = parameters.simulation('S_slope');
    S_threshold     = parameters.simulation('S_threshold');
    S_to_C_alpha    = parameters.simulation('S_to_C_alpha'); % learning rate
    S_sigma         = V_sigma; % (deg) receptive field size

    % C
    C_N             = S_N*R_N;
    C_activation    = zeros(1,C_N);
    C_firingrate    = zeros(1,C_N);
    
    C_tau           = parameters.simulation('C_tau'); % (s)
    C_to_R_alpha    = parameters.simulation('C_to_R_alpha');
    C_to_R_psi      = parameters.simulation('C_to_R_psi');
    C_w_INHB        = 0/C_N;
    C_slope         = parameters.simulation('C_slope');
    C_threshold     = parameters.simulation('C_threshold');
    
    % Allocate buffer space
    V_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    R_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    S_firing_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
    C_firing_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);

    V_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    R_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    S_activation_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
    C_activation_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);
    
    %% Simulate
    totalTicID = tic;
    for epoch=1:numEpochs,
        
        if isTraining,
            disp(['Starting epoch #' num2str(epoch)]);
        end
        
        for period=1:numPeriods
            
            periodTicID = tic;
            
            % Get number of time steps in this period
            numTimeSteps = timeStepsInPeriod(period);
            
            maxNumberOfVisibleTargets = length(stimuli.stimuli{period}.headCenteredTargetLocations);
    
            % Load Stimuli for period
            retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces;
            saccadeTimes        = stimuli.stimuli{period}.saccadeTimes;
            saccadeTargets      = stimuli.stimuli{period}.saccadeTargets;
            numSaccades         = stimuli.stimuli{period}.numSaccades;

            % Setup static working variables
            R_preference_comparison_matrix = repmat(R_preferences, maxNumberOfVisibleTargets, 1); % used to compute driving term in V
            
            if(numSaccades > 0),
                saccade_times = repmat(saccadeTimes,S_N,1); % Used to get F
                saccade_time_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                saccade_target_offset = exp(-((repmat(S_preferences', 1, numSaccades) - repmat(saccadeTargets, S_N, 1)).^2)./(2*S_sigma^2)); % term multiplied by F
            end
            
            % Reset network variables
            periodSaveCounter = 1;
            
            % do later?
    
            % Run period
            for t=1:numTimeSteps,

                % Turn time step into real time
                time = stepToTime(t);

                % Activation

                % R
                R_inhibition = R_w_INHB*sum(R_firingrate);
                C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate');
                V_to_R_excitation = V_to_R_psi*V;
                R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation' - R_inhibition + V_to_R_excitation);

                % V
                retinalTargets = retinalTargetTraces(:,t);
                diff = R_preference_comparison_matrix - repmat(retinalTargets,1,R_N);
                diff(isnan(diff)) = inf; % for nan values, make inf, so that exponantiation gives no contribution
                gauss = exp((-(diff).^2)./(2*V_sigma^2));

                V_driver = V_psi*sum(gauss,1); % add upp all targets for each neuron to get one driving V value
                V = V + (dt/V_tau)*(-V + V_driver);

                % C
                C_inhibition = C_w_INHB*sum(C_firingrate);
                %R_to_C_excitation = R_to_C_psi*(R_to_C_weights*R_firingrate');
                V_to_C_excitation = V_to_C_psi*(V_to_C_weights*V');
                S_to_C_excitation = S_to_C_psi*(S_to_C_weights*S_firingrate');
                C_activation = C_activation + (dt/C_tau)*(-C_activation + V_to_C_excitation' + S_to_C_excitation' - C_inhibition); % _to_C_excitation

                % S
                if(numSaccades > 0),
                    F = (saccade_time_offset <= time) & (time <= saccade_times); % check both conditions: y-z <= x <= y
                    S_driver = S_psi*sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
                    S_activation = S_activation + (dt/S_tau)*(-S_activation + S_driver');
                else
                    S_activation = S_activation + (dt/S_tau)*(-S_activation); 
                end
                

                % Weight Update
                if isTraining,

                    C_to_R_weights = C_to_R_weights + dt*C_to_R_alpha*(R_firingrate'*C_firingrate);
                    S_to_C_weights = S_to_C_weights + dt*S_to_C_alpha*(C_firingrate'*S_firingrate);
                    R_to_C_weights = R_to_C_weights + dt*R_to_C_alpha*(C_firingrate'*R_firingrate);

                    % Normalize
                    C_to_R_norm = 1./sqrt(squeeze(sum(C_to_R_weights.^2))); 
                    C_to_R_weights = bsxfun(@times,C_to_R_weights,C_to_R_norm);

                    S_to_C_norm = 1./sqrt(squeeze(sum(S_to_C_weights.^2))); 
                    S_to_C_weights = bsxfun(@times,S_to_C_weights,S_to_C_norm); 

                    R_to_C_norm = 1./sqrt(squeeze(sum(R_to_C_weights.^2))); 
                    R_to_C_weights = bsxfun(@times,R_to_C_weights,R_to_C_norm);
                end

                % Compute firing rates
                R_firingrate = 1./(1 + exp(-2*R_slope*(R_activation - R_threshold)));
                S_firingrate = 1./(1 + exp(-2*S_slope*(S_activation - S_threshold)));
                C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_threshold)));

                % Save activity                
                if (~isTraining || isTraining && parameters.saveActivityInTraining) && mod(t, outputSavingRate) == 0,
                    
                    V_firing_history(:, periodSaveCounter, period, epoch) = V;
                    R_firing_history(:, periodSaveCounter, period, epoch) = R_firingrate;
                    S_firing_history(:, periodSaveCounter, period, epoch) = S_firingrate;
                    C_firing_history(:, periodSaveCounter, period, epoch) = C_firingrate;

                    V_activation_history(:, periodSaveCounter, period, epoch) = V;
                    R_activation_history(:, periodSaveCounter, period, epoch) = R_activation;
                    S_activation_history(:, periodSaveCounter, period, epoch) = S_activation;
                    C_activation_history(:, periodSaveCounter, period, epoch) = C_activation;
                    
                    
                    % Count one more dt
                    periodSaveCounter = periodSaveCounter + 1;
                end 
            end
            
            periodFinishTime = toc(periodTicID);
            disp(['Finished period #' num2str(period) ' of ' num2str(numPeriods) ' in ' num2str(periodFinishTime) 's.']);
            
        end
        
        % Output intermediate network
        if mod(epoch, parameters.saveNetworksAtEpochMultiples) == 0,
            
            disp('Saving trained network to disk...');
            save([simulationFolder filesep 'TrainedNetwork_e' num2str(epoch) '.mat'] , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights');
        end
    end
    
    totalFinishTime = toc(totalTicID);
    disp(['Completed ' num2str(numEpochs) ' epochs in ' num2str(totalFinishTime) 's.']);
    
    % Output final network
    if isTraining,
        disp('Saving trained network to disk...');
        save([simulationFolder filesep 'TrainedNetwork.mat'] , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights', 'R_N', 'S_N', 'C_N');
    end
    
    % Output activity
    disp('Saving activity...');
    save([simulationFolder filesep 'activity' outputpostfix '.mat'] , 'V_firing_history' ...
                                                                    , 'R_firing_history' ...
                                                                    , 'S_firing_history' ...
                                                                    , 'C_firing_history' ...
                                                                    , 'V_activation_history' ...
                                                                    , 'R_activation_history' ...
                                                                    , 'S_activation_history' ...
                                                                    , 'C_activation_history' ...
                                                                    , 'numEpochs' ...
                                                                    , 'numPeriods' ...
                                                                    , 'R_N', 'S_N', 'C_N');
    
    disp('Done...');
    
    function r = stepToTime(i)
        r = (i-1)*dt;
    end

end