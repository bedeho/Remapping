
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
            [StimuliControl_Neurons, StimuliControl_indexes] = AnalyzeStimuliControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'StimuliControl_Neurons', 'StimuliControl_indexes');
            
        elseif strcmp(type,'SaccadeControl'),
            
            disp('Doing saccade control task analysis...');
            [saccade_response] = AnalyzeSaccadeControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'saccade_response');
            
        elseif strcmp(type,'DuhamelRemapping'),
            
            
            disp('Doing duhamel remapping task analysis...');
            [DuhamelRemapping_Neurons, DuhamelRemapping_indexes] = AnalyzeDuhamelRemapping(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemapping_Neurons', 'DuhamelRemapping_indexes');
            
        elseif strcmp(type,'DuhamelRemappingTrace'),
            
            
            disp('Doing duhamel remapping trace task analysis...');
            [DuhamelRemappingTrace_Neurons, DuhamelRemappingTrace_indexes] = AnalyzeDuhamelRemappingTrace(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemapping_Neurons', 'DuhamelRemapping_indexes');
            
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
    latency = [StimuliControl_Neurons(:).latency];
    maxLatency = max(latency);
    minLatency = min(latency);
    dh = (minLatency - minLatency)/11;
    
    if(dh == 0),
        
        lat = latency(1);
        
        left_bins = fliplr(lat:-0.01:(lat-0.05));
        right_bins = lat:0.01:(lat+0.05);
        
        x = [left_bins right_bins(2:end)]
    else
        x = minLatency:dh:maxLatency;
    end
    
    bar(x,hist(latency,x));
    
    xlabel('Time (s)');
    ylabel('Frequency');

    saveas(f,[netDir filesep 'StimuliControl-summary.png']);
    close(f);
    
    %% Duhamel remapping analysis
    % stim: [StimuliControl_Neurons, StimuliControl_indexes]
    % duhamel: [DuhamelRemapping_analyzedNeurons, DuhamelRemapping_indexes]
    
    f = figure;
    hold on;
    
    % Iteratea all neurons studied in remapping context
    for index=DuhamelRemapping_indexes,
    
        % Find the stimuli control latency of this neuron
        j = find(StimuliControl_indexes == index);
        
        if(length(j) == 1)
            plot(StimuliControl_Neurons(j).latency, DuhamelRemapping_Neurons(index).latency,'ro');
        else
            error('STIM CONTROL TASK SHOULD ONLY TEST EACH NEURON ONCE.');
        end
    end
    
    
    
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
    
    xlabel('Stimulus Control Latency (s)');
    ylabel('Remapping Latency (s)');
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    plot([-0.5 0.5],[-0.5 0.5],'--b');
    axis square
    
    saveas(f,[netDir filesep 'DuhamelRemapping-summary.png']);
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
    
    %% BACKUP
    
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
end