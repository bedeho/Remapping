
%
%  Remapping.m
%  Remapping
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Remapping(sourcefolder, outputFolder, paramters,enablePlasticity)

    % Dynamical quantities
    Duration = 2; % (s)
    dt = 0.010; % (s)
    numTimeSteps = ceil(Duration/dt);

    % Random seed
    rng(77);

    %% Stimuli
    
    % Visual
    headCenteredTargetLocations = [15]; % (deg)
    maxNumberOfVisibleTargets = length(headCenteredTargetLocations);
    
    targetOffIntervals{1} = []%[0.6 2];%[0.5 0.9;]; % (s) [start_OFF end_OFF; start_OFF end_OFF] <==== Make dt multiples
    targetOffIntervals{2} = [];%[0.1 0.2;];
    
    assert(length(targetOffIntervals) >= maxNumberOfVisibleTargets, 'On off history not provided for all targets.');
    
    % Saccadic
    initialEyePosition = 0; % (deg), init = 13 deg
    saccadeSpeed = 300; % (deg/s)
    saccadeTimes = [1.2]; % (s) % <==== Make dt multiples, ALLOW FOR SUFFIEIENCT INTERSACCADE TIME TO COMPLETE SACCADES!
    saccadeTargets = [20]; % (deg)
    numSaccades = length(saccadeTimes);
    
    assert(length(saccadeTimes) >= length(saccadeTargets), 'Number of saccade times and targets must match.');

    % Allocate space for traces
    % Index i => Time (i-1)*dt
    eyePositionTrace = zeros(1, numTimeSteps);
    eyePositionTrace(1) = initialEyePosition;

    % Generate eye position trace
    numCompletedSaccades = 0;
    for t=2:numTimeSteps,
        
        presentTime = stepToTime(t);
        
        % Have we completed all saccades beginning before time timestep t?
        if numCompletedSaccades == numSaccades || presentTime <= saccadeTimes(numCompletedSaccades+1),

            % yes, so lets just continue fixating
            eyePositionTrace(t) = eyePositionTrace(t-1);

        else
            % no, there was a not completed saccade beginning before
            % timestep t
            
            % was this the first step across this saccade onset time?
            previousTimeSteptime = stepToTime(t-1);
            timeOffset = abs(previousTimeSteptime - saccadeTimes(numCompletedSaccades+1));
            
            if(timeOffset < dt),
                
                % yes it was, so mini saccade starting
                % point needs to correct for delay between prior time step
                % and saccade onset event to find saccade start position
                
                reduceSaccadeTimeInPresentTimeStepWith = timeOffset;
            else
                
                % no, we passed at some previous time, so mini saccade starting
                % point is just last eye position
                
                reduceSaccadeTimeInPresentTimeStepWith = 0;
            end
            
            % Will one more stime step saccading from
            % miniSaccadeStartPosition take us past saccade target?
            
            saccadeOffset = eyePositionTrace(t-1) - saccadeTargets(numCompletedSaccades+1);
            timeStepSaccadeMagnitude = (dt - reduceSaccadeTimeInPresentTimeStepWith)*saccadeSpeed;
            
            if (timeStepSaccadeMagnitude > saccadeOffset),
                
                % yes, so lets not do the whole thing
                eyePositionTrace(t) = saccadeTargets(numCompletedSaccades+1);
                
                % and lets saccade as completed
                numCompletedSaccades = numCompletedSaccades + 1;
                
            else
                
                % no, then lets do the full thing and keep going.
                eyePositionTrace(t) = eyePositionTrace(t-1) + -1*sign(saccadeOffset)*saccadeSpeed*dt;
                
            end
            
        end
        
    end

    % Generate retinal target location trace
    retinalTargetTraces = zeros(maxNumberOfVisibleTargets, numTimeSteps);
    for h=1:maxNumberOfVisibleTargets,
        
        % Make trace: r = h - e
        retinalTargetTraces(h, :) = headCenteredTargetLocations(h) - eyePositionTrace;
        
        % Cancel out parts where target is not present
        offIntervals = targetOffIntervals{h};
        [numIntervals x] = size(offIntervals);
        for i=1:numIntervals,
            
            % Get interval
            interval = offIntervals(i,:);
            
            % Translate from time to timesteps
            timeStepInterval = timeToTimeStep(interval);
            
            % Cancel out, i.e. not visible
            retinalTargetTraces(h, timeStepInterval(1):timeStepInterval(2)) = nan;
        end
    end
    
    function r = stepToTime(i)
        r = (i-1)*dt;
    end
    
    function i = timeToTimeStep(t)
        i = floor(t/dt) + 1;
    end

    %% Dynamics
    
    % LIP: Retina
    R_eccentricity = 45;
    R_preferences = -R_eccentricity:1:R_eccentricity;
    R_N = length(R_preferences);
    
    R_tau = 0.100; % (s)
    C_to_R_psi = 0.08; % 0.15
    R_w_INHB = 0; %0.7
    
    V_tau = 0.400; % (s)
    V_psi = 4;
    V_sigma = 5; % (deg) receptive field size
    V = zeros(1,R_N);
    
    C_to_R_alpha = 0.1; % learning rate
    
    R_slope = 1;
    R_threshold = 2.0;
    
    R_activation = zeros(1,R_N);
    R_firingrate = zeros(1,R_N);
    R_history = zeros(R_N,numTimeSteps);
    R_history2 = zeros(R_N,numTimeSteps);
    
    % FEF: Saccade Plan
    S_eccentricity = 30;
    S_preferences = -S_eccentricity:1:S_eccentricity;
    S_N = length(S_preferences);
    S_delay_sigma = 0.4; % (s)
    S_presaccadicOffset = S_delay_sigma*randn(1,S_N);
    S_presaccadicOffset(S_presaccadicOffset < 0) = -S_presaccadicOffset(S_presaccadicOffset < 0);
    
    S_tau = 0.300; % (s)
    S_psi = 1;
    S_sigma = V_sigma; % (deg) receptive field size
    
    S_slope = 6;
    S_threshold = 0.2;
    
    S_activation = zeros(1,S_N);
    S_firingrate = zeros(1,S_N);
    S_history = zeros(S_N,numTimeSteps);
    S_history2 = zeros(S_N,numTimeSteps);
    
    % SC?: Comb
    C_N = S_N*R_N;
    C_tau = 0.100; % (s)
    R_to_C_psi = 1.0; % 4
    S_to_C_psi = 6.2;
    C_w_INHB = 0/C_N;
    
    R_to_C_alpha = 0.1; % learning rate
    S_to_C_alpha = 0.1; % learning rate
    
    C_slope = 50;
    C_threshold = 1.0;
    
    C_activation = zeros(1,C_N);
    C_firingrate = zeros(1,C_N);
    C_history = zeros(C_N,numTimeSteps);
    C_history2 = zeros(C_N,numTimeSteps);
    
    %% Setup Weights
    
    % Random
    C_to_R_weights = rand(R_N,C_N);
    R_to_C_weights = rand(C_N,R_N);
    S_to_C_weights = rand(C_N,S_N);
    
    %{
    % Hardwire
    [X Y Z] = meshgrid(R_preferences, S_preferences, R_preferences);
    
    % C_to_R_weights
    C_to_R_raw = exp(-(((Z - (X-Y))).^2)./(2*V_sigma^2));
    C_to_R_reshaped = reshape(C_to_R_raw,C_N,R_N)';
    C_to_R_norm = 1./sqrt(squeeze(sum(C_to_R_reshaped.^2))); 
    C_to_R_weights = bsxfun(@times,C_to_R_reshaped,C_to_R_norm); % Normalize
    
    % R_to_C_weights
    R_to_C_raw = exp(-((Z - X).^2)./(2*V_sigma^2));
    R_to_C_reshaped = reshape(R_to_C_raw,C_N,R_N);
    R_to_C_norm = 1./sqrt(squeeze(sum(R_to_C_reshaped.^2))); 
    R_to_C_weights = bsxfun(@times,R_to_C_reshaped,R_to_C_norm); % Normalize
    
    % S_to_C_weights
    [X Y Z] = meshgrid(R_preferences, S_preferences, S_preferences);
    
    S_to_C_raw = exp(-((Z - Y).^2))./(2*V_sigma^2);
    S_to_C_reshaped = reshape(S_to_C_raw, C_N,S_N);
    S_to_C_norm = 1./sqrt(squeeze(sum(S_to_C_reshaped.^2))); 
    S_to_C_weights = bsxfun(@times,S_to_C_reshaped,S_to_C_norm); % Normalize
    %}
    
    %% Integrate
    
    % Forward euler: dY = Y + (dt/tau)*(f(Y,t))
    %
    % Computation order respecting dependancy
    % 1. Activation = old (activations, firing rates, weights
    % 2. Weight = old firing rates
    % 3. Firing rate = old activations
    
    % static working variables
    R_preference_comparison_matrix = repmat(R_preferences, maxNumberOfVisibleTargets, 1); % used to compute driving term in V
    
    % meshgrid alternative: even tighter
    %[saccade_times saccade_presaccadicoffset] = meshgrid(saccadeTimes, S_presaccadicOffset);
    %saccade_time_offset = saccade_times - saccade_presaccadicoffset;
    
    saccade_times = repmat(saccadeTimes,S_N,1); % Used to get F
    saccade_time_offset = saccade_times - repmat(S_presaccadicOffset',1,numSaccades); % Used to get F
    saccade_target_offset = exp(-((repmat(S_preferences', 1, numSaccades) - repmat(saccadeTargets, S_N, 1)).^2)./(2*S_sigma^2)); % term multiplied by F
    
    tic
    for t=1:numTimeSteps,
        
        % Turn time step into real time
        time = stepToTime(t);
        
        %% Activation
        
        % R
        R_inhibition = R_w_INHB*sum(R_firingrate);
        C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate'); 
        R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation' - R_inhibition + V);
        
        % V
        retinalTargets = retinalTargetTraces(:,t);
        diff = R_preference_comparison_matrix - repmat(retinalTargets,1,R_N);
        diff(isnan(diff)) = inf; % for nan values, make inf, so that exponantiation gives no contribution
        gauss = exp((-(diff).^2)./(2*V_sigma^2));
        
        V_driver = V_psi*sum(gauss,1); % add upp all targets for each neuron to get one driving V value
        V = V + (dt/V_tau)*(-V + V_driver);
        
        % C
        C_inhibition = C_w_INHB*sum(C_firingrate);
        R_to_C_excitation = R_to_C_psi*(R_to_C_weights*R_firingrate');
        S_to_C_excitation = S_to_C_psi*(S_to_C_weights*S_firingrate');
        C_activation = C_activation + (dt/C_tau)*(-C_activation + R_to_C_excitation' + S_to_C_excitation' - C_inhibition);
        
        % S
        F = (saccade_time_offset <= time) & (time <= saccade_times); % check both conditions: y-z <= x <= y
        S_driver = S_psi*sum(bsxfun(@times, F, saccade_target_offset),2); % cannot be done with matrix mult since exponential depends on i
        S_activation = S_activation + (dt/S_tau)*(-S_activation + S_driver'); 
        
        %% Weight Update
        if enablePlasticity,
            
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
        
        %% Firing rates & Save history
        R_firingrate = 1./(1 + exp(-2*R_slope*(R_activation - R_threshold)));
        S_firingrate = 1./(1 + exp(-2*S_slope*(S_activation - S_threshold)));
        C_firingrate = 1./(1 + exp(-2*C_slope*(C_activation - C_threshold)));
        
        % Save in history
        R_history(:,t) = R_firingrate; % R_firingrate, R_activation
        S_history(:,t) = S_firingrate;
        C_history(:,t) = C_firingrate; %C_firingrate; C_activation
        
        R_history2(:,t) = R_activation; % R_firingrate, R_activation
        S_history2(:,t) = S_activation;
        C_history2(:,t) = C_activation;
        
    end
    
    toc
    
    % Plot
    figure('Position', [100, 100, 1049, 895]);

    subplot(4,2,1);
    imagesc(flipud(R_history));
    colorbar
    
    subplot(4,2,2);
    imagesc(flipud(R_history2));
    colorbar
    
    subplot(4,2,3);
    imagesc(flipud(S_history));
    colorbar
    
    subplot(4,2,4);
    imagesc(flipud(S_history2));
    colorbar
    
    subplot(4,2,5);
    imagesc(flipud(C_history));
    colorbar
    
    subplot(4,2,6);
    imagesc(flipud(C_history2));
    colorbar
    
    subplot(4,2,7);
    plot(eyePositionTrace, 'r');
    hold on;
    plot(retinalTargetTraces' , 'b');
    xlabel('Time step');
    legend({'Eye Position','Stimuli Retinal Locations'})
    %ylim(min(min(eyePositionTrace),min(min(retinalTargetTraces)) max(max(eyePositionTrace),max(max(retinalTargetTraces)))]);
end