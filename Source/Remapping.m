
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
    C_to_R_weights_dilutionmap = network.C_to_R_weights_dilutionmap;
    S_to_C_weights_dilutionmap = network.S_to_C_weights_dilutionmap;
    V_to_C_weights_dilutionmap = network.V_to_C_weights_dilutionmap;

    %% Load dynamical parameters
    
    % R  =======================================
    R_preferences   = parameters.simulation('R_preferences');
    R_N             = size(C_to_R_weights,1); %length(R_preferences);
    R_activation    = zeros(1,R_N);
    R_firingrate    = zeros(1,R_N);
    
    R_tau           = parameters.simulation('R_tau');
    R_w_INHB        = parameters.simulation('R_w_INHB');
    R_slope         = parameters.simulation('R_slope');
    R_threshold     = parameters.simulation('R_threshold');

    R_covariance_threshold = parameters.simulation('R_covariance_threshold');

    % K  =======================================
    K_tau           = parameters.simulation('K_tau');
    P_tau           = parameters.simulation('P_tau');
    P_psi           = parameters.simulation('P_psi');
    
    K_I_psi         = parameters.simulation('K_I_psi');
    K_onset_delays     = parameters.simulation('K_onset_delays');
    K_supression_delay = parameters.simulation('K_supression_delay');
    
    %{
    if(maxNumberOfVisibleTargets > 0),
        K = zeros(maxNumberOfVisibleTargets, R_N);
    else
        K = 0;
    end
    %}
    
    K                = zeros(1, R_N);
    P                = zeros(1, R_N);
    
    K_delaysTimeSteps = durationToSteps(K_onset_delays, dt);
    K_supression_delayTimeSteps = durationToSteps(K_supression_delay, dt);
    
    % V =======================================
    V_sigma         = parameters.simulation('V_sigma');
    V_to_C_psi      = parameters.simulation('V_to_C_psi');
    V_to_C_alpha    = parameters.simulation('V_to_C_alpha'); % learning rate
    V               = zeros(1,R_N);
    V_firingrate    = zeros(1,R_N);
    
    V_supression_delay = parameters.simulation('V_supression_delay');
    V_supression_delayTimeSteps = durationToSteps(V_supression_delay, dt);
    
    % S  =======================================
    S_preferences   = parameters.simulation('S_preferences');
    S_N             = size(S_to_C_weights, 2);%length(S_preferences);
    
    S_activation    = zeros(1,S_N);
    S_firingrate    = zeros(1,S_N);
    
    S_presaccadicOffset = parameters.simulation('S_presaccadicOffset'); %
    S_tau           = parameters.simulation('S_tau'); % (s)
    S_to_C_psi      = parameters.simulation('S_to_C_psi');
    S_psi           = parameters.simulation('S_psi');
    S_to_C_alpha    = parameters.simulation('S_to_C_alpha'); % learning rate
    S_sigma         = V_sigma; % (deg) receptive field size
    
    S_presaccadic_onset = parameters.simulation('S_presaccadic_onset')*ones(1,S_N);
    S_trace_length      = parameters.simulation('S_trace_length')*ones(1,S_N);

    % C  =======================================
    C_N             = size(S_to_C_weights, 1); % S_N*R_N;
    C_activation    = zeros(1,C_N);
    C_firingrate    = zeros(1,C_N);
    
    C_tau           = parameters.simulation('C_tau'); % (s)
    C_to_R_alpha    = parameters.simulation('C_to_R_alpha');
    C_to_R_psi      = parameters.simulation('C_to_R_psi');
    C_w_INHB        = parameters.simulation('C_w_INHB');
    C_slope         = parameters.simulation('C_slope');
    C_threshold     = parameters.simulation('C_threshold');
    
    %C_thresholds    = parameters.simulation('C_thresholds');
    
    % K thresholds
    %g = randn(1, C_N);
    %g(g > 0) = 0; % high limit: mean
    %g(g < -1) = -1; % low limit: mean - one std
    %C_thresholds   = parameters.simulation('C_threshold_sigma')*g + parameters.simulation('C_threshold');
            
    
    % Allocate buffer space
    saveOutput = (~isTraining) || (isTraining && parameters.saveActivityInTraining);
    
    if(saveOutput),
        
        V_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
        R_firing_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
        S_firing_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
        C_firing_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);

        V_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
        R_activation_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
        S_activation_history = zeros(S_N, numSavedTimeSteps, numPeriods, numEpochs);
        C_activation_history = zeros(C_N, numSavedTimeSteps, numPeriods, numEpochs);

        extra_history = zeros(R_N, numSavedTimeSteps, numPeriods, numEpochs);
    end
    
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
            Duration = stepToTime(numTimeSteps,dt);
            
            maxNumberOfVisibleTargets = length(stimuli.stimuli{period}.headCenteredTargetLocations);
    
            % Load Stimuli for period
            retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces;
            saccadeTimes        = stimuli.stimuli{period}.saccadeTimes;
            saccadeTargets      = stimuli.stimuli{period}.saccadeTargets;
            stimOnsetTimes      = stimuli.stimuli{period}.stimOnsetTimes;
            
            stimOnsetTimeSteps  = timeToTimeStep(stimOnsetTimes,dt);
            saccOnsetTimeSteps  = timeToTimeStep(saccadeTimes,dt);

            numSaccades         = length(stimuli.stimuli{period}.saccadeTimes);
            numStimOnsetTimes   = length(stimOnsetTimes);

            % Setup for indicator function
            targetOffIntervals = stimuli.stimuli{period}.targetOffIntervals{1};
            
            tmp = targetOffIntervals;
            tmp = tmp';
            tmp = tmp(:)';
            
            if(~isempty(stimOnsetTimes)),
                tmp = [0 tmp Duration];
                tmp = reshape(tmp, 2, numel(tmp)/2);
                on_off_timestep_pairs = timeToTimeStep(tmp',dt);
            end
            
            if(~isempty(targetOffIntervals)),
            
                % offset in P
                if(targetOffIntervals(1,1) == 0),
                    off_time_steps = timeToTimeStep(targetOffIntervals(2:end,1), dt);
                else
                    off_time_steps = timeToTimeStep(targetOffIntervals(:,1), dt);
                end
            
            else
                off_time_steps = [];
            end
            
            
            % Setup static working variables
            R_preference_comparison_matrix = repmat(R_preferences, maxNumberOfVisibleTargets, 1);

            stimOnset_comparison_matrix = repmat(stimOnsetTimeSteps', 1, R_N); 
            
            K_delay_comparison_matrix = repmat(K_delaysTimeSteps, numStimOnsetTimes, 1);
            
            %stimulus_presence_indicator = ~isnan(retinalTargetTraces);
            
            if(numSaccades > 0),
                saccade_times = repmat(saccadeTimes,S_N,1); % Used to get F
                saccade_time_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                
                % classic
                %saccade_time_pos_offset = saccade_times + repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                %saccade_time_neg_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
                
                % assymmetry
                saccade_time_pos_offset = saccade_times + repmat(S_trace_length',1,numSaccades);
                saccade_time_neg_offset = saccade_times - repmat(S_presaccadic_onset',1,numSaccades);
                
                saccade_target_offset = exp(-((repmat(S_preferences', 1, numSaccades) - repmat(saccadeTargets, S_N, 1)).^2)./(2*S_sigma^2)); % term multiplied by F
            end
            
            % Reset network variables
            periodSaveCounter   = 2; % First dt is automatically saved
            V                   = V*0;
            K                   = K*0; % must be reallocated since teh number of visiible targets can in theory change, although it never does in practice.
            P                   = P*0;
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
                V_old = V;
                K_old = K;
                P_old = P;
                
                % Prepare response of all neurons to all visible targets
                if(~isempty(retinalTargetTraces)),
                    
                    % retinal locations of stimuli in the precedng time step
                    retinalTargets = retinalTargetTraces(:,t-1);
                    flat_gauss = visul_response(retinalTargets);

                    % clipping of retinal trace
                    indicator = zeros(1,R_N);
                    for q=1:R_N,
                        
                        starts = on_off_timestep_pairs(:,1) + K_delaysTimeSteps(q);
                        stops = on_off_timestep_pairs(:,2);
                        in_interval = (starts <= (t-1)) & ((t-1) <= stops);
                        indicator(q) = any(in_interval);
                    end                  

                    K_feed_forward = indicator.*flat_gauss;
                else
                    flat_gauss = 0;
                    K_feed_forward = 0;
                end
                
                % Visual onset spike
                if(~isempty(stimOnsetTimes)),
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    
                    visual_onset = zeros(1, R_N);
                    visual_onset(yes_delta_event) = flat_gauss(yes_delta_event);
                else
                    visual_onset = 0;
                end
                
                % R ======================================= 
                
                % P
                if(~isempty(stimOnsetTimeSteps) && any(precedingTimeStep==off_time_steps)),
                    offset_decay = P_psi*K_old;
                else
                    offset_decay = 0;
                end
                
                if(~isempty(saccOnsetTimeSteps) && any(precedingTimeStep==saccOnsetTimeSteps+K_supression_delayTimeSteps)),
                    P_sacc_supression = P_old; %* K_supress
                    K_sacc_supression = K_old;
                    
                    %R_sacc_supression = R_activation;
                    
                    offset_decay = 0; %phenomenological trick in case both happen at the same time, could rewrite system, but why bother?
                else
                    P_sacc_supression = 0;
                    K_sacc_supression = 0;
                    
                    %R_sacc_supression = 0;
                end
                
                %P_visual_input = K_I_psi*K_feed_forward;
                P = P_old + (dt/P_tau)*(-P_old) - P_sacc_supression + offset_decay; % P_visual_input
                
                % K
                K_visual_onset = K_I_psi*visual_onset;
                K_visual_input = K_I_psi*K_feed_forward;
                K = K_old + (dt/K_tau)*(-K_old + K_visual_input + P_old) + K_visual_onset - K_sacc_supression; %
                
                % R
                %R_visual_onset = K_psi*visual_onset;
                R_global_inhibition = R_w_INHB*sum(R_firingrate); % R_activation.*
                C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate')';
                R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation - R_global_inhibition + K ); % + R_visual_onset, - R_background

                % V =======================================
                
                % stim
                if(~isempty(stimOnsetTimeSteps)),
                    
                    if(any(precedingTimeStep==stimOnsetTimeSteps))
                        V_onset = flat_gauss;
                    else
                        V_onset = 0;
                    end
                    
                    % sacc onset/supression
                    if(~isempty(saccOnsetTimeSteps) && any(precedingTimeStep==saccOnsetTimeSteps+V_supression_delayTimeSteps) && (stimOnsetTimeSteps < saccOnsetTimeSteps)),
                        V_sacc_supression = V_old;
                        V_sacc_excitation = flat_gauss;
                    else
                        V_sacc_supression = 0;
                        V_sacc_excitation = 0;
                    end
                    
                else
                    V_onset = 0;
                    
                    V_sacc_supression = 0;
                    V_sacc_excitation = 0;
                end
                
                V = V_old + V_onset - V_sacc_supression + V_sacc_excitation; % (dt/V_tau)*(-V_old) +

                % C =======================================
                
                C_inhibition = C_w_INHB*sum(C_firingrate);
                V_to_C_excitation = V_to_C_psi*(V_to_C_weights*V_firingrate');
                S_to_C_excitation = S_to_C_psi*(S_to_C_weights*S_firingrate');
                C_activation = C_activation + (dt/C_tau)*(-C_activation + V_to_C_excitation' + S_to_C_excitation' - C_inhibition);

                % S =======================================
                if(numSaccades > 0),
                    F = (saccade_time_neg_offset <= time_precedingdt) & (time_precedingdt <= saccade_time_pos_offset); % check both conditions: y-z <= x <= y
                    %F = (saccade_time_offset <= time_precedingdt) & (time_precedingdt <= saccade_times); % check both conditions: y-z <= x <= y
                    %S_driver = S_psi*sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
                    S_driver = S_psi*sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
                    S_activation = S_activation + (dt/S_tau)*(-S_activation + S_driver');
                else
                    S_activation = S_activation + (dt/S_tau)*(-S_activation); 
                end
                
                %% Weight Update
                if isTraining,

                    % Learning rule
                    %C_to_R_weights = C_to_R_weights + dt*C_to_R_alpha*(1.5
                    %- C_to_R_weights).*(R_firingrate'*C_firingrate); % Bounded learning 
                    
                    % CLASSIC
                    %C_to_R_weights = C_to_R_weights + dt*C_to_R_alpha*(R_firingrate'*C_firingrate);
                    
                    % COVARIANCE
                    C_to_R_weights = C_to_R_weights + dt*C_to_R_alpha*((R_firingrate' - R_covariance_threshold)*C_firingrate);
                    C_to_R_weights(C_to_R_weights < 0) = 0;
                    
                    S_to_C_weights = S_to_C_weights + dt*S_to_C_alpha*(C_firingrate'*S_firingrate);
                    V_to_C_weights = V_to_C_weights + dt*V_to_C_alpha*(C_firingrate'*V_firingrate);

                    % Diluted connectivity
                    C_to_R_weights = C_to_R_weights.*C_to_R_weights_dilutionmap;
                    S_to_C_weights = S_to_C_weights.*S_to_C_weights_dilutionmap;
                    V_to_C_weights = V_to_C_weights.*V_to_C_weights_dilutionmap;

                    %{
                    figure;
                    subplot(1,2,1);
                    imagesc(C_to_R_weights);

                    % Normalize
                    C_to_R_weights = normalizeWeightVector(C_to_R_weights);
                    S_to_C_weights = normalizeWeightVector(S_to_C_weights);
                    V_to_C_weights = normalizeWeightVector(V_to_C_weights);

                    subplot(1,2,2);
                    imagesc(C_to_R_weights);
                    %}
                    
                    % Normalize
                    C_to_R_weights = normalizeWeightVector(C_to_R_weights);
                    S_to_C_weights = normalizeWeightVector(S_to_C_weights);
                    V_to_C_weights = normalizeWeightVector(V_to_C_weights);
                    
                end

                %% Compute firing rates
                R_firingrate = 1./(1 + exp(-2*R_slope*(R_activation - R_threshold)));
                %R_firingrate = R_activation;
                
                S_firingrate = S_activation;
                
                C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_threshold)));
                
                %C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_thresholds)));
                
                %C_dynamic_threshold = prctile(C_activation, C_percentile);

                %C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_dynamic_threshold - C_threshold)));
                
                %V_firingrate = 1./(1 + exp(-2*C_slope*(V - V_threshold))); 
                V_firingrate = V; 
                
                %% Save activity                
                if (saveOutput),
                    
                    V_firing_history(:, periodSaveCounter, period, epoch) = V_firingrate;
                    R_firing_history(:, periodSaveCounter, period, epoch) = R_firingrate;
                    S_firing_history(:, periodSaveCounter, period, epoch) = S_firingrate; % S_activation;
                    C_firing_history(:, periodSaveCounter, period, epoch) = C_firingrate;

                    V_activation_history(:, periodSaveCounter, period, epoch) = V;
                    R_activation_history(:, periodSaveCounter, period, epoch) = R_activation;
                    S_activation_history(:, periodSaveCounter, period, epoch) = S_activation;
                    C_activation_history(:, periodSaveCounter, period, epoch) = C_activation;
                    
                    extra_history(:, periodSaveCounter, period, epoch) = C_to_R_excitation;%ones(size(P))*C_dynamic_threshold;
                    
                    % Count one more dt
                    periodSaveCounter = periodSaveCounter + 1;
                end 
                
                %{
                % Visualize
                if(strcmp(stimuliName,'basic-Training_Coordinated') && epoch >= 0 && t >= 40 ), % && t==141 period == 6, period == numPeriods, && 

                    
                    n=130;

                    %figure;
                    %plot(C_to_R_weights(40, :));
                    %axis tight
                    
                    %figure;
                    
                    subplot(6,3,1);
                    plot(S_to_C_weights(n,:));
                    axis tight;
                    ylim([0 0.5]);
                    title('S->C');

                    subplot(6,3,2);
                    plot(V_to_C_weights(n,:));
                    axis tight;
                    ylim([0 0.5]);
                    title('V->C');

                    subplot(6,3,3);
                    plot(C_to_R_weights(:, n));
                    axis tight;
                    %ylim([0 0.4]);
                    title('C->R');

                    subplot(6,3,[4 5 6]);
                    imagesc(S_firing_history(:, :, period, epoch));
                    hold on;
                    if ~isempty(saccOnsetTimeSteps), plot([saccOnsetTimeSteps saccOnsetTimeSteps],[ones(S_N,1) S_N*ones(S_N,1)],'r'); end
                    title('S');

                    subplot(6,3,[7 8 9]);
                    imagesc(V_firing_history(:, :, period, epoch));
                    hold on;
                    if ~isempty(saccOnsetTimeSteps), plot([saccOnsetTimeSteps saccOnsetTimeSteps],[ones(R_N,1) S_N*ones(R_N,1)],'r'); end
                    title('V');

                    subplot(6,3,[10 11 12]);
                    imagesc(R_firing_history(:, :, period, epoch));
                    hold on;
                    if ~isempty(saccOnsetTimeSteps), plot([saccOnsetTimeSteps saccOnsetTimeSteps],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
                    title('R');

                    subplot(6,3,[13 14 15]);
                    cla
                    plot(C_firing_history(n, :, period, epoch));
                    hold on;
                    if ~isempty(saccOnsetTimeSteps), plot([saccOnsetTimeSteps saccOnsetTimeSteps],[0 1.1],'r'); end
                    title('C Firing of neuron n');

                    t_show = 50;%60;

                    subplot(6,3,[16 17 18]);
                    cla
                    tmp = C_activation_history(:, t_show, period, epoch);
                    [nelements, centers] =hist(tmp(:),100);
                    bar(centers,nelements);
                    hold on;
                    n_activation = C_activation_history(n, t_show, period, epoch);
                    plot([n_activation n_activation],[0 max(nelements)],'r','LineWidth',2);
                    plot([C_threshold C_threshold],[0 max(nelements)],'g','LineWidth',2);
                    
                    axis tight;
                    title('C activation');

                end
                %}
                
            end

            % HACK to save C for C Probe Task, we dont have space to do
            % this properly
            if numSaccades > 0 && ~isTraining,

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
        % 'R_to_R_excitatory_weights', 'R_to_R_inhibitory_weights',
        save([simulationFolder filesep 'TrainedNetwork.mat'] , 'C_to_R_weights', 'S_to_C_weights', 'V_to_C_weights',  'C_to_R_weights_dilutionmap', 'S_to_C_weights_dilutionmap', 'V_to_C_weights_dilutionmap', 'R_N', 'S_N', 'C_N'); %'V_to_R_weights'
    end
    
    if saveOutput && numel(C_firing_history) > 2000000/2, % dont save if bigger than 100MB
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
                                                                            , 'V_firing_history' ...
                                                                            , 'R_firing_history' ...
                                                                            , 'S_firing_history' ...
                                                                            , 'C_firing_history' ...
                                                                            , 'V_activation_history' ...
                                                                            , 'R_activation_history' ...
                                                                            , 'S_activation_history' ...
                                                                            , 'C_activation_history' ...
                                                                            , 'C_firing_history_flat' ...
                                                                            , 'C_activation_history_flat' ...
                                                                            , 'extra_history');
                                                                        
    disp('Done...');
              
    function flat_gauss = visul_response(ret_target_locations)

        diff = R_preference_comparison_matrix - repmat(ret_target_locations,1,R_N);
        diff(isnan(diff)) = inf; % for nan values, make inf, so that exponantiation gives no contribution
        gauss = exp((-(diff).^2)./(2*V_sigma^2)); % gauss has dimensions: maxNumberOfVisibleTargets X R_N
        flat_gauss = sum(gauss,1); % collapse stimulus dimension, so that each neuron has one driving sum of exponentials, one exponential per stimulus

    end
            
end