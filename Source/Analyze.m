
%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function analysisSummary = Analyze(netDir, stimulinames)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;

    for i = 1:length(stimulinames),
        
        % Load stimuli
        disp(['Loading stimuli: ' stimulinames{i}]);
        stimuliFile = [STIMULI_FOLDER stimulinames{i} filesep 'stim.mat'];
        stimuli  = load(stimuliFile);
        type = stimuli.stimulitype;
        
        % Load activity
        disp(['Loading activity: ' stimulinames{i}]);
        activityFile = [netDir filesep 'activity-' stimulinames{i} '.mat'];
        activity = load(activityFile);
        %activity = LoadActivity(activityFile);
        
        if strcmp(type,'StimuliControl'),
            
            disp('Doing stimuli control task analysis...');
            [StimuliControl_Result] = AnalyzeStimuliControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'StimuliControl_Result');
            
        elseif strcmp(type,'SaccadeControl'),
            
            disp('Doing saccade control task analysis...');
            [SaccadeControl_Result] = AnalyzeSaccadeControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'SaccadeControl_Result');
            
        elseif strcmp(type,'DuhamelRemapping'),
            
            disp('Doing duhamel remapping task analysis...');
            [DuhamelRemapping_Result] = AnalyzeDuhamelRemapping(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemapping_Result');
            
        elseif strcmp(type,'DuhamelRemappingTrace'),
            
            disp('Doing duhamel remapping trace task analysis...');
            [DuhamelRemappingTrace_Result] = AnalyzeDuhamelRemapping(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemappingTrace_Result');
            
        elseif strcmp(type,'DuhamelTruncation'),
            
            [DuhamelTruncation_Neurons, DuhamelTruncation_indexes] = AnalyzeDuhamelTruncation(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelTruncation_Neurons', 'DuhamelTruncation_indexes');
              
        elseif strcmp(type,'Kusonoki'),
            
            disp('Doing Kusonoki analysis...');
            [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activity, stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'kusonokiSTIMAlignedAnalysis', 'kusonokiSACCAlignedAnalysis');
            
        elseif strcmp(type,'CLayerProbe'),
            
            disp('Doing C Layer Probe analysis...');
            [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity,  stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'CLabeProbe_Neurons_V', 'CLabeProbe_Neurons_V');
            
        else
            disp(['Unkonwn stimuli: ' num2str(stimulinames{i})]);
        end
            
    end
    
    %% Stimuli Control
    f = figure;
    latency = [StimuliControl_Result(:).latency];
    maxLatency = max(latency);
    minLatency = min(latency);
    dh = (minLatency - minLatency)/11;
    
    if(dh == 0),
        
        lat = latency(1);
        
        left_bins = fliplr(lat:-0.01:(lat-0.05));
        right_bins = lat:0.01:(lat+0.05);
        
        x = [left_bins right_bins(2:end)];
    else
        x = minLatency:dh:maxLatency;
    end
    
    bar(x,hist(latency,x));
    
    xlabel('Time (s)');
    ylabel('Frequency');

    saveas(f,[netDir filesep 'StimuliControl-summary.png']);
    close(f);
    
    %% Duhamel remapping analysis
    
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
    
    %% Start plotting
    
    % 1. scatter remap latency vs. stim control latency
    f = figure;
    hold on;
    plot(remapping_latency, stim_latency, 'or');
    plot([-0.5 0.5],[-0.5 0.5],'--b'); % y=x bar
    
    xlabel('Stimulus Control Latency (s)');
    ylabel('Remapping Latency (s)');
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    axis square;
    
    saveas(f,[netDir filesep 'DuhamelRemapping-summary.png']);
    close(f);

    % 2. scatter stim index. vs sacc index.
    f = figure;
    hold on;
    plot(sacc_index, stim_index, 'or');
    plot([0 0],[-1 1],'--g'); % x=0 bar
    plot([-1 1],[0 0],'--g'); % y=0 bar
    
    xlabel('Saccade Index');
    ylabel('Stimulus Index');
    xlim([-1 1]);
    ylim([-1 1]);
    axis square;
    
    saveas(f,[netDir filesep 'DuhamelRemapping-summary-2.png']);
    close(f);
    
    % 3. remapping index distibution
    f = figure;
    hold on;
    x = 0:0.1:sqrt(2);
    bar(x,hist(remapping_index,x));
    
    xlabel('Remapping Index');
    ylabel('Frequency');
    axis square;
    
    saveas(f,[netDir filesep 'DuhamelRemapping-summary-3.png']);
    close(f);
    
    %% Duhamel trace remapping analysis
    
    %% Duhamel truncation analysis
    
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
    
    %% C PRobe
    f = figure;
    plot(CLabeProbe_Neurons_V, CLabeProbe_Neurons_S, 'or');
    xlabel('Retinal Locaton');
    ylabel('Saccade Location');
    
    saveas(f,[netDir filesep 'CLayerProbe-summary.png']);
    
    close(f);
    
    analysisSummary = 0;
    
    function remappingAnalysis(remapping_indexes, remapping_neurons, stimcontrol_neurons, sacccontrol_neurons)
        

    end

end