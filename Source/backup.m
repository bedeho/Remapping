    %{
    % Load Duhamel Truncation file
    duhamelTruncationFile = [experimentFolder filesep 'analysis-basic-DuhamelTruncation.mat'];
    if(exist(duhamelTruncationFile, 'file')),
        
        duhamelTruncation = load(duhamelTruncationFile);
        DuhamelRemappingTrace_Result = duhamelRemapping.DuhamelRemappingTrace_Result;
        
        %{
        f = figure;
        hold on;
        plot([DuhamelTruncation_Result(:).saccadeonset_response], [DuhamelTruncation_Result(:).stim_stim_offset_response], 'or');
        plot([0 1],[0 1],'--b'); % y=x bar
        xlabel('Truncation Saccade Onset Response');
        ylabel('Stimulus Offset Response');
        xlim([0 1]);
        ylim([0 1]);
        axis square;

        saveas(f,[netDir filesep 'DuhamelTruncation-summary.png']);
        close(f);
        %}
    end
    %}


%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

    %{
    %% Stimuli Control
    
    f = figure;
    latency = [StimuliControl_Result(:).latency];
    maxLatency = max(latency);
    minLatency = min(latency);
    dh = (maxLatency - minLatency)/21;
    
    if(dh == 0),
        
        lat = latency(1);
        
        left_bins = fliplr(lat:-0.01:(lat-0.05));
        right_bins = lat:0.01:(lat+0.05);
        
        x = [left_bins right_bins(2:end)];
    else
        x = minLatency:dh:maxLatency;
    end
    
    % test in case experiment fails, analysis crashed
    if ~isnan(x),
        bar(x,hist(latency,x));
    end
    
    xlabel('Time (s)');
    ylabel('Frequency');

    saveas(f,[netDir filesep 'StimuliControl-summary.png']);
    close(f);
    
    %{
    %}
    
    %% Duhamel remapping analysis plotting
    remappingAnalysis(DuhamelRemapping_Result, 'DuhamelRemapping');
    
    %% Duhamel trace remapping analysis
    remappingAnalysis(DuhamelRemappingTrace_Result, 'DuhamelRemappingTrace');
    
    %% Duhamel truncation analysis
    f = figure;
    hold on;
    plot([DuhamelTruncation_Result(:).saccadeonset_response], [DuhamelTruncation_Result(:).stim_stim_offset_response], 'or');
    plot([0 1],[0 1],'--b'); % y=x bar
    xlabel('Truncation Saccade Onset Response');
    ylabel('Stimulus Offset Response');
    xlim([0 1]);
    ylim([0 1]);
    axis square;

    saveas(f,[netDir filesep 'DuhamelTruncation-summary.png']);
    close(f);
    
    
    %% Kusonoki
    if(exist('kusonokiSTIMAlignedAnalysis') && exist('kusonokiSACCAlignedAnalysis')),
        
        f = figure;

        arr_stim = kusonokiSTIMAlignedAnalysis;
        subplot(2,1,1);
        hold on;
        errorbar([arr_stim(:).current_mean], [arr_stim(:).current_std],'-or');
        errorbar([arr_stim(:).future_mean], [arr_stim(:).future_std],'-ob');
        legend('Current RF Trials','Future RF Trials');
        ylim([0 1]);

        arr_sacc = kusonokiSACCAlignedAnalysis;
        subplot(2,1,2);
        hold on;
        errorbar([arr_sacc(:).current_mean], [arr_sacc(:).current_std],'-or');
        errorbar([arr_sacc(:).future_mean], [arr_sacc(:).future_std],'-ob');
        legend('Current RF Trials','Future RF Trials');
        ylim([0 1]);

        saveas(f,[netDir filesep 'Kusonoki-summary.png']);
        close(f);
    end
    
    
    %% C PRobe
    f = figure;
    plot(CLabeProbe_Neurons_V, CLabeProbe_Neurons_S, 'or');
    xlabel('Retinal Locaton');
    ylabel('Saccade Location');
    
    saveas(f,[netDir filesep 'CLayerProbe-summary.png']);
    
    close(f);
    
    
    function remappingAnalysis(remapping_result, name)
        
        % 1. scatter remap latency vs. stim control latency
        f = figure;
        hold on;
        plot([remapping_result(:).stimLatency], [remapping_result(:).remappingLatency], 'or');
        plot([-0.2 0.2],[-0.2 0.2],'--b'); % y=x bar

        xlabel('Stimulus Control Latency (s)');
        ylabel('Remapping Latency (s)');
        xlim([-0.2 0.2]);
        ylim([-0.2 0.2]);
        axis square;

        saveas(f,[netDir filesep name '-summary.png']);
        close(f);

        % 2. scatter stim index. vs sacc index.
        f = figure;
        hold on;
        plot([remapping_result(:).sacc_index], [remapping_result(:).stim_index], 'or');
        plot([0 0],[-1 1],'--g'); % x=0 bar
        plot([-1 1],[0 0],'--g'); % y=0 bar

        xlabel('Saccade Index');
        ylabel('Stimulus Index');
        xlim([-1 1]);
        ylim([-1 1]);
        axis square;

        saveas(f,[netDir filesep name '-summary-0.png']);
        close(f);


        % 3. remapping index distibution
        f = figure;
        hold on;
        x = 0:0.1:sqrt(2);
        bar(x,hist([remapping_result(:).remapping_index],x));

        xlabel('Remapping Index');
        ylabel('Frequency');
        axis square;

        saveas(f,[netDir filesep name '-summary-2.png']);
        close(f);
        
    end
    %}
    
    
    
    
    
    
    
    
    
    
    

%{
                % CLASSIC - SOM based
                if(~isempty(stimOnsetTimes)),
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    
                    visual_onset = zeros(1, R_N);
                    visual_onset(yes_delta_event) = K_psi*flat_gauss(yes_delta_event);
                else
                    visual_onset = 0;
                end
                
                K_visual_input = K_I_psi*K_feed_forward;
                K = K_old + (dt/K_tau)*(-K_old + K_visual_input);
                
                
                C_to_R_excitation       = C_to_R_psi*(C_to_R_weights*C_firingrate')';
                
                R_SOM_inhibition        = R_neg_attractor_psi*(R_to_R_inhibitory_weights*R_firingrate')';
                
                R_global_inhibition     = R_w_INHB*sum(R_firingrate);
                
                R_attractor             = R_attractor_psi*(R_to_R_excitatory_weights*R_firingrate')';
                
                R_activation            = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation - R_global_inhibition + K_visual_input + R_attractor + R_SOM_inhibition) + visual_onset;
                %}


            %{
            on_off_time_pairs = zeros(numStimOnsetTimes, 2);
            for j=1:numStimOnsetTimes,
                
                % on time step
                on_off_timestep_pairs(j, 1) = timeToTimeStep(stimOnsetTimes(j), dt);
                
                % off time step
                offset_pair = stimuli.stimuli{period}.targetOffIntervals{1}; % we use the same on off for all stim! I dont even give an *ff at this point
                
                if(~isempty(offset_pair)),
                    
                    
                    
                    if(offset_pair(1,1) > 0), % stimuli was immediately turned off
                        
                        % if so
                    
                    else
                        
                    end
                    
                    
                    on_off_timestep_pairs(j, 2) = timeToTimeStep(offset_pair(1,2), dt);
                else
                    on_off_timestep_pairs(j, 2) = numTimeSteps;
                end
            end
            %}


%
%  GenerateExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.


                
                %{
                R_SOM_inhibition        = R_neg_attractor_psi*((R_to_R_inhibitory_weights*R_activation')');
                R_attractor             = R_attractor_psi*(R_to_R_weights*R_activation')';
                C_to_R_excitation       = C_to_R_psi*(C_to_R_weights*C_firingrate')';
                R_visual                = R_psi*flat_gauss;
                
                R_net_input = C_to_R_excitation + R_SOM_inhibition + R_visual + R_attractor;
                R_sig = 1./(1 + exp(-2*R_slope*(R_net_input - R_threshold)));
                
                R_activation = R_activation + (dt/R_tau)*(-R_activation + R_sig);%; - R_activation*K_sacc_supression;
                %}

                %R_inhibition_scaling = R_activation;
                %R_inhibition_scaling(R_inhibition_scaling > 1) = 1; % do linear rectified function of R scaling.


                %C_to_R_inhibition = 0;%R_inhibition_scaling*C_to_R_psi_neg*sum(C_firingrate);

                %
                %R_SOM_inhibition = R_inhibition_scaling.*(R_neg_attractor_psi*(R_to_R_inhibitory_weights*R_firingrate')');
                
                %R_attractor = R_attractor_psi*(R_to_R_excitatory_weights*R_firingrate')';


                %K_tau_dynamic = R_tau_rise + (flat_gauss <= R_tau_threshold)*(R_tau_decay-R_tau_rise);
                %K = K_old + (dt./K_tau_dynamic).*(-K_old + R_psi*flat_gauss) + K_onset_spike - K_sacc_supression;
                
                %K = K_old + (dt/K_tau)*(-K_old) + K_onset_spike - K_old*K_sacc_supression;
                
                % V =======================================
                
                % E & 
                %E_tau_dynamic = E_tau_rise + exp(-(flat_gauss.^2)/(2*E_tau_sigma^2))*(E_tau_decay-E_tau_rise);
                %E = E_old + (dt./E_tau_dynamic).*(-E_old + flat_gauss);
                %V = V_old + (dt/V_tau)*(-V_old + E_to_V_psi*E_old) - V_old*K_sacc_supression;
                
                
                
                
                
                
                
                
                
                
                
                

%{

%
%  GenerateExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

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
%}

%
%  Remapping.m
%  Remapping
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.

%{
K_gauss = flat_gauss;

K_tau_dynamic = K_tau*ones(1,R_N);%R_tau_rise + (K_gauss <= R_tau_threshold)*(R_tau_decay-R_tau_rise);

if(~isempty(stimOnsetTimes)),

    % Time comparison for delta impulse must be done in
    % time steps, not continous time, otherwise we it will
    % almost surely miss delta(0)
    precedingTimeStep = t-1;

    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
    delta_sum = sum(delta == 0, 1);
    yes_delta_event = (delta_sum > 0);
    no_delta_event = ~yes_delta_event;

    % Update neurons with delta event
    K(yes_delta_event) = K_old(yes_delta_event) + K_psi*K_gauss(yes_delta_event); % dirac delta case, +K_psi*delta_sum(yes_delta_event).*flat_gauss(yes_delta_event)

    % Update neurons without delta event: standard FE
    K(no_delta_event) = K_old(no_delta_event) + (dt./K_tau_dynamic(no_delta_event)).*(-K(no_delta_event)); % + K_gauss(no_delta_event), continous case when there is no dirac delta event
else

    % Do all neuron at ones
    K = K_old + (dt./K_tau_dynamic).*(-K_old); % K_gauss
end
%}                


%{
                K_gauss = flat_gauss;
                
                K_tau_dynamic = K_tau*ones(1,R_N);%R_tau_rise + (K_gauss <= R_tau_threshold)*(R_tau_decay-R_tau_rise);
                
                if(~isempty(stimOnsetTimes)),
                    
                    % Time comparison for delta impulse must be done in
                    % time steps, not continous time, otherwise we it will
                    % almost surely miss delta(0)
                    precedingTimeStep = t-1;
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    no_delta_event = ~yes_delta_event;

                    % Update neurons with delta event
                    K(yes_delta_event) = K_old(yes_delta_event) + K_psi*K_gauss(yes_delta_event); % dirac delta case, +K_psi*delta_sum(yes_delta_event).*flat_gauss(yes_delta_event)

                    % Update neurons without delta event: standard FE
                    K(no_delta_event) = K_old(no_delta_event) + (dt./K_tau_dynamic(no_delta_event)).*(-K(no_delta_event)); % + K_gauss(no_delta_event), continous case when there is no dirac delta event
                else
                    
                    % Do all neuron at ones
                    K = K_old + (dt./K_tau_dynamic).*(-K_old); % K_gauss
                end
%}

%{
                K_gauss = flat_gauss;
                
                K_tau_dynamic = K_tau*ones(1,R_N);%R_tau_rise + (K_gauss <= R_tau_threshold)*(R_tau_decay-R_tau_rise);
                
                if(~isempty(stimOnsetTimes)),
                    
                    % Time comparison for delta impulse must be done in
                    % time steps, not continous time, otherwise we it will
                    % almost surely miss delta(0)
                    precedingTimeStep = t-1;
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    no_delta_event = ~yes_delta_event;

                    % Update neurons with delta event
                    K(yes_delta_event) = K_old(yes_delta_event) + K_psi*K_gauss(yes_delta_event); % dirac delta case, +K_psi*delta_sum(yes_delta_event).*flat_gauss(yes_delta_event)

                    % Update neurons without delta event: standard FE
                    K(no_delta_event) = K_old(no_delta_event) + (dt./K_tau_dynamic(no_delta_event)).*(-K(no_delta_event)); % + K_gauss(no_delta_event), continous case when there is no dirac delta event
                else
                    
                    % Do all neuron at ones
                    K = K_old + (dt./K_tau_dynamic).*(-K_old); % K_gauss
                end
                
%}                
                


  
                %{
                
                R_inhibition = R_w_INHB*sum(R_firingrate);
                C_to_R_excitation = C_to_R_psi*(C_to_R_weights*C_firingrate');
                %V_to_R_excitation = V_to_R_psi*V;
                
                % R asymmetric decay
                %R_visual_excitation = R_psi*gauss; % classic
                R_visual_excitation = R_psi*flat_gauss; %sum(gauss,1);
                R_total_exication = C_to_R_excitation' + R_visual_excitation;

                R_total_input = R_total_exication - R_inhibition;

                % classic dynamic tau
                %R_tau_dynamic = R_tau_rise + exp(-(R_total_exication.^2)/(2*R_tau_sigma^2))*(R_tau_decay-R_tau_rise);

                %R_tau_dynamic = R_tau_rise + (R_total_exication <= R_tau_threshold)*(R_tau_decay-R_tau_rise);

                R_tau_dynamic = R_tau_rise + (abs(R_total_input) <= R_tau_threshold)*(R_tau_decay-R_tau_rise);

                if(~isempty(stimOnsetTimes)),
                    
                    % Time comparison for delta impulse must be done in
                    % time steps, not continous time, otherwise we it will
                    % almost surely miss delta(0)
                    precedingTimeStep = t-1;
                    
                    delta = (precedingTimeStep - stimOnset_comparison_matrix - K_delay_comparison_matrix);
                    delta_sum = sum(delta == 0, 1);
                    yes_delta_event = (delta_sum > 0);
                    no_delta_event = ~yes_delta_event;

                    % Update neurons with delta event
                    R_activation(yes_delta_event) = R_activation(yes_delta_event) + R_psi*K_psi*delta_sum(yes_delta_event).*R_visual_excitation(yes_delta_event);

                    % Update neurons without delta event: standard FE
                    R_activation(no_delta_event) = R_activation(no_delta_event) + (dt./R_tau_dynamic(no_delta_event)).*(-R_activation(no_delta_event) + R_total_input(no_delta_event));
                else
                    
                    % Do all neuron at ones
                    R_activation = R_activation + (dt./R_tau_dynamic).*(-R_activation + R_total_exication - R_inhibition);
                end
                
                %}
                
                %{
                CLASSIC: This is with variable onset times for different
                stimuli

                if(~isempty(stimOnsetTimes)),
                    
                    % Time comparison for delta impusle must be done in
                    % time steps, not continous time, otherwise we it will
                    % almost surely miss delta(0)
                    precedingTimeStep = t-1;
                    
                    delta = precedingTimeStep - stimOnset_comparison_matrix + K_delay_comparison_matrix;
                    delta_sum = sum(delta == 0, 1);

                    % We have to iterate neurons since delta function needs
                    % explicit formula
                    for i=1:R_N,

                        if(delta_sum(i) > 0),
                            K(i) = K_old(i) + K_psi*delta_sum(i);
                        else
                            K(i) = K_old(i) + (dt/K_tau)*(-K_old(i) + stimulus_presence_indicator(:,t-1));
                        end
                    end

                else
                    K = 0;
                end
                
                R_visual_excitation = sum(K.*gauss,1);
                R_total_exication = C_to_R_excitation' + R_visual_excitation;
                R_tau = R_tau_rise + exp(-(R_total_exication.^2)/(2*R_tau_sigma^2))*(R_tau_decay-R_tau_rise);
                R_activation = R_activation + (dt./R_tau).*(-R_activation + R_total_exication - R_inhibition );
                %}

                % gauss drives R
                %R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation' - R_inhibition + R_psi*gauss);
                
                % E drives R
                %R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation' - R_inhibition + E_to_R_psi*E);
                
                %classic: 
                %R_activation = R_activation + (dt/R_tau)*(-R_activation + C_to_R_excitation' - R_inhibition + V_to_R_excitation);
                
            
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% backup

    %% Baseline response
    
    %{
    
    % Get time steps in question
    baselineTimeSteps   = 1:(onsetTimeStep-1);
    
    % Extract the given time steps from all neurons in all periods
    baselineActivity    = R_firing_history(:, baselineTimeSteps, :, 1);
    
    % Integrate to find response
    baselineResponse    = squeeze(trapz(baselineActivity,2)); % Integrate
    
    % Should be all identical across neurons and periods, so we pick the first period for all of them.
    baselineResponse    = baselineResponse(:,1);
    
    % Normaliztion step, gives normalized (sp/s) units to response
    baselineResponse    = baselineResponse/(length(baselineTimeSteps)-1);
    
    %% Stimulus response
    
    % Get time steps in question
    activityTimeSteps   = timeToTimeStep(stimuliOnsetDelay + responseWindowStart:dt:responseWindowEnd, dt);
    
    % Extract the given time steps from all neurons in all periods
    stim_activity       = R_firing_history(:, activityTimeSteps, :, 1); % [onsetTimeStep+50:250]
    
    % Integrate to find response
    stim_response       = squeeze(trapz(stim_activity,2));
    
    % Normaliztion step, gives normalized (sp/s) units to response
    stim_response       = stim_response/(length(activityTimeSteps) - 1);
    
    %% Receptive field location
    
    % Level normalization
    stim_response_modenormalized = stim_response - repmat(mode(stim_response, 2),1, numPeriods);
    stim_response_modenormalized(stim_response_modenormalized < 0) = 0;
    
    % Center of mass
    normalization_response = sum(stim_response_modenormalized,2);
    location = (stim_response*stimuli.headCenteredTargetLocations')./normalization_response;
    
    %}
    
    
    %{
    
        activityTimeSteps   = timeToTimeStep(stimuliOnsetDelay + (responseWindowStart:dt:responseWindowEnd), dt);  % Get time steps in question
        stim_activity       = neuronActivity(activityTimeSteps); % [onsetTimeStep+50:250]
        stim_response       = squeeze(trapz(stim_activity,2)); % Integrate to find response
        stim_response       = stim_response/(length(activityTimeSteps) - 1); % Normaliztion step, gives normalized (sp/s) units to response
        
                baselineTimeSteps = 1:(onsetTimeStep-1);
        baselineActivity = neuronActivity(baselineTimeSteps);
        baselineResponse = squeeze(trapz(baselineActivity,2)); % Integrate
        baselineResponse = baselineResponse/(length(baselineTimeSteps)-1); % Normaliztion step, gives normalized (sp/s) units to response
        
                offsetTimeSteps     = timeToTimeStep(stimuliOffsetPeriod + stimuliOnsetDelay + (0:dt:responseWindowDuration), dt);  % Get time steps in question
        offset_activity     = neuronActivity(offsetTimeSteps); % [onsetTimeStep+50:250]
        offset_response     = squeeze(trapz(offset_activity,2)); % Integrate to find response
        offset_response     = offset_response/(length(offsetTimeSteps) - 1); % Normaliztion step, gives normalized (sp/s) units to response
    %}
    
    
    
    
 %{   
%
%  AnalyzeDuhamelRemappingTrace.m
%  Remapping
%
%  Created by Bedeho Mender on 10/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [DuhamelRemappingTrace_Neurons, DuhamelRemappingTrace_indexes] = AnalyzeDuhamelRemappingTrace(activity, stimuli)

    error('No longer in use.');



    % Check if this is manual run 
    if nargin == 0,
        
        disp('Loading input files...');
        %LoadActivity
        activity = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-DuhamelRemappingTrace.mat');
        stimuli  = load('/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Stimuli/basic-DuhamelRemappingTrace/stim.mat');
    end
    
    % Get data
    R_firing_history = activity.R_firing_history;
    
    % Set parameters
    dt                   = activity.dt;
    R_eccentricity       = stimuli.R_eccentricity;
    numPeriods           = activity.numPeriods;
    numEpochs            = activity.numEpochs;
    saccadeOnset         = stimuli.saccadeOnset;
    FutureRFLocations    = stimuli.FutureRFLocations;
    
    % Analysis params
    %latencyWindowSize   = 0.020; % (s), colby papers
    %latencyWindowLength = ceil(latencyWindowSize/dt);
    RF_inclusion_th     = 5; % (deg) neurons this far away from any given trial are analysed together
    responseWindowDuration = 0.200;
    %responseThreshold   = 0.5;
    
    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    
    c = 1;
    DuhamelRemapping_neuronIndexes = [];
    
    for p=1:numPeriods,
        
        % Target location
        futureRFLocation = FutureRFLocations(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,futureRFLocation - RF_inclusion_th):1:min(R_eccentricity,futureRFLocation + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data neuron
            neuronActivity  = R_firing_history(neuronIndex, :, p, 1);
            
            % Offset response
            saccadeonset_response = normalizedIntegration(neuronActivity, dt, saccadeOnset, responseWindowDuration);
            
            %% DEBUG - looks GOOD!
            %figure;plot(neuronActivity);hold on; plot(timeToTimeStep([saccadeOnset saccadeOnset], dt), [0 1], 'r');
            
            % Save
            DuhamelRemappingTrace_Neurons(c).index                   = neuronIndex;
            DuhamelRemappingTrace_Neurons(c).saccadeonset_response   = saccadeonset_response;
            
            DuhamelRemappingTrace_indexes(c) = neuronIndex;
            
            c = c + 1;
            
        end
    end


    
end
    %}
    
    
        
    %{
    old style, with a single neuron included multiple times because
    futureRF is close to , but not identical to its location
        c = 1;
    DuhamelRemapping_neuronIndexes = [];
    
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation_HeadCentered = stimuli.headCenteredTargetLocations(p);
        futureRFLocation = currentRFLocation_HeadCentered - stimuli.saccadeTargets(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,futureRFLocation - RF_inclusion_th):1:min(R_eccentricity,futureRFLocation + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data for best period of each neuron
            responseVector  = R_firing_history(neuronIndex, :, p, 1);
            
            % Find latency and duration
            [latencyTimeStep, duration] = findNeuronalLatency(responseThreshold, responseVector, latencyWindowLength);
            
            %figure;plot(responseVector);hold on; plot([latencyTimeStep latencyTimeStep],[0 1], 'r');
            
            % Save
            DuhamelRemapping_Neurons(c).index            = neuronIndex;
            DuhamelRemapping_Neurons(c).latency          = stepToTime(latencyTimeStep, dt)-saccadeOnset;
            DuhamelRemapping_Neurons(c).Duration         = duration*dt;
            
            DuhamelRemapping_indexes(c) = neuronIndex;
            
            c = c + 1;
            
        end
    end
    %}
    
    
    %% BACKUP
    
    %
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

    
    % This was in the stimuli analysis case: what does it do?
    
                %{
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , ...
                    'baselineResponse', ...
                    'stim_response', ...
                    'location', ...
                    'foundOnset', ...
                    'foundOffset', ...
                    'latencyTimeStep', ...
                    'durationTimeStep');
                            
            R_N = size(latencyTimeStep,2);                
            f = figure;
            imagesc(neuronResponse);
            ylabel('Neuron');
            xlabel('Time');
            hold on;
            plot(latencyTimeStep,1:R_N,'wo');
            saveas(f,[netDir filesep stimulinames{i} '.png']);
            close(f);
            %}
    
    
        
    %{
    old, when neurons were included multiple times,
    but it picks neurons that were only included ones,
    rather than pick the best instances of all neurons ever
    included, so update this if it is ever used again s
    for i=1:length(StimuliControl_indexes),
        
        index_1 = StimuliControl_indexes(i);
        j = find( == index_1);
        
        if(length(j) == 1),
            
            plot(StimuliControl_Neurons(i).latency, DuhamelRemapping_Neurons(j).latency,'ro');
            %disp('found duhamel remapping neuron');
        end
        
    end
    %}
    
%
%  GenerateExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%
    
        %{
    % Num neurons
    numNeurons = length(DuhamelRemapping_Result);
    
    StimuliControl_rf      = [StimuliControl_Result(:).receptiveField];
    SaccadeControl_saccade = [SaccadeControl_Result(:).receptiveField];
    
    % Iteratea all neurons studied in remapping context
    stim_index   = zeros(1, numNeurons);
    sacc_index   = zeros(1, numNeurons);
    stim_latency = zeros(1, numNeurons);

    for i=1:numNeurons,

        % future RF neuron index
        futureRf = DuhamelRemapping_Result(i).futureRF;
        
        % saccade executed
        saccade = DuhamelRemapping_Result(i).saccade;
        
        % original location prior to saccade
        currentRf = futureRf + saccade;
        
        % Find the stimuli control latency of this neuron
        j_stim = find(StimuliControl_rf == currentRf);
        j_sacc = find(SaccadeControl_saccade == saccade);

        % Check that we only get one hit
        if(length(j_stim) ~= 1 || length(j_sacc) ~= 1)
            error('STIM & SACC. CONTROL TASK SHOULD ONLY TEST EACH NEURON ONCE, AND MUST TEST ALL NEURONS.');
        end
        
        % Compute stimulus index: based on stim onset response in control
        % task and remapping task
        stim_index(i) = DuhamelRemapping_Result(i).saccadeonset_response - StimuliControl_Result(j_stim).stimulus_response;
        
        % Compute saccade index: based on saccade onset response in control
        % task and remapping task
        sacc_index(i) = DuhamelRemapping_Result(i).saccadeonset_response - SaccadeControl_Result(j_sacc).saccadeonset_response;
        
        % Get latency
        stim_latency(i) = StimuliControl_Result(j_stim).latency;

    end

    remapping_index     = sqrt(stim_index.^2 + sacc_index.^2); % according to L.M.Heiser,Colby (2006)
    remapping_latency   = [DuhamelRemapping_Result(:).latency];
    DuhamelRemapping_Index = [DuhamelRemapping_Result(:).index];
    %}


%
%  AnalyzeDuhamelTruncation.m
%  Remapping
%
%  Created by Bedeho Mender on 10/07/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%


%{

    assert(numEpochs == 1, 'There is more than one epoch, hence this is not a testing stimuli');
    
    %% Latency & Duration
    
    c = 1;
    DuhamelTruncation_indexes = [];
    
    for p=1:numPeriods,
        
        % Target location
        currentRFLocation_HeadCentered = stimuli.headCenteredTargetLocations(p);
        
        % Find neurons that are close enough
        neuron_RFLocations = max(-R_eccentricity,currentRFLocation_HeadCentered - RF_inclusion_th):1:min(R_eccentricity,currentRFLocation_HeadCentered + RF_inclusion_th);
        
        for f=neuron_RFLocations,
        
            % Find neuron
            neuronIndex = R_eccentricity + f + 1;

            % Get data neuron
            neuronActivity  = R_firing_history(neuronIndex, :, p, 1);
            
            % Offset response
            saccadeonset_response = normalizedIntegration(neuronActivity, dt, saccadeOnset, responseWindowDuration);
            
            %% DEBUG - looked good
            %figure;plot(neuronActivity);hold on; plot(timeToTimeStep([saccadeOnset saccadeOnset], dt), [0 1], 'r');
            
            % Save
            DuhamelTruncation_Neurons(c).index                   = neuronIndex;
            DuhamelTruncation_Neurons(c).saccadeonset_response   = saccadeonset_response;
            
            DuhamelTruncation_indexes(c) = neuronIndex;
            
            c = c + 1;
            
        end
    end
%}

%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%


    %% Kusonoki
    
    % Get a plot up and running? MANUAL HOME BREW ANALYTICS
    %{

    %activity = load(activityFile);
    activity = LoadActivity(activityFile);
    stimuli  = load(stimuliFile);

    period = 20;
    R_N = activity.R_N;
    dt = activity.dt;

    f=figure;
    x = activity.R_firing_history(:, :, period, 1);

    s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);

    nr = stimuli.stimuli{period}.stimOnsetNr;
    stim = timeToTimeStep(stimuli.stimulusOnsetTimes(nr), dt);
    stim_off = timeToTimeStep(stimuli.stimulusOnsetTimes(nr) + stimuli.stimulusDuration, dt);

    imagesc(x);

    hold on;

    plot([stim stim],[ones(R_N,1) R_N*ones(R_N,1)],'--w','LineWidth', 4); % STIM ON

    plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'w','LineWidth', 4); % saccade

    plot([stim_off stim_off],[ones(R_N,1) R_N*ones(R_N,1)],'--w','LineWidth', 4); % STIM OFF

    title(num2str(period));

    saveas(f,[netDir filesep 'summary.png']);
    close(f);
    %}