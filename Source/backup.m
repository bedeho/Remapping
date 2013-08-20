
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