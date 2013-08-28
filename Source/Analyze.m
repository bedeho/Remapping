
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
    
    % keep track of thesee
    stim_control_activity  = [];
    sacc_control_activity = [];
    stim_stimuli = [];
    sacc_stimuli = [];
    
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
            StimuliControl_Result = AnalyzeStimuliControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'StimuliControl_Result');
            
            stim_stimuli = stimuli;
            stim_control_activity = activity;
            stim_control_results = StimuliControl_Result;
            
        elseif strcmp(type,'SaccadeControl'),
            
            disp('Doing saccade control task analysis...');
            SaccadeControl_Result = AnalyzeSaccadeControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'SaccadeControl_Result');
            
            sacc_stimuli = stimuli;
            sacc_control_activity = activity;
            
        elseif strcmp(type,'DuhamelRemapping'),
            
            disp('Doing duhamel remapping task analysis...');
            DuhamelRemapping_Result = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemapping_Result');
            
        elseif strcmp(type,'DuhamelRemappingTrace'),
            
            disp('Doing duhamel remapping trace task analysis...');
            DuhamelRemappingTrace_Result = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemappingTrace_Result');
            
        elseif strcmp(type,'DuhamelTruncation'),
            
            DuhamelTruncation_Result = AnalyzeDuhamelTruncation(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelTruncation_Result');
              
        elseif strcmp(type,'Kusonoki'),
            
            disp('Doing Kusonoki analysis...');
            [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activity, stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'kusonokiSTIMAlignedAnalysis', 'kusonokiSACCAlignedAnalysis');
            
        elseif strcmp(type,'CLayerProbe'),
            
            disp('Doing C Layer Probe analysis...');
            [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity,  stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'CLabeProbe_Neurons_S', 'CLabeProbe_Neurons_V');
            
        else
            disp(['Unkonwn stimuli: ' num2str(stimulinames{i})]);
        end
            
    end
    
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
    
    bar(x,hist(latency,x));
    
    xlabel('Time (s)');
    ylabel('Frequency');

    saveas(f,[netDir filesep 'StimuliControl-summary.png']);
    close(f);
    
    
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
    f = figure;
    hold on;
    
    arr = kusonokiSTIMAlignedAnalysis;
    %arr = kusonokiSACCAlignedAnalysis;
    
    errorbar([arr(:).current_mean], [arr(:).current_std],'-or');
    errorbar([arr(:).future_mean], [arr(:).future_std],'-ob');
    
    legend('Current RF Trials','Future RF Trials');
    ylim([0 0.6]);

    saveas(f,[netDir filesep 'Kusonoki-summary.png']);
    close(f);

    %% C PRobe
    f = figure;
    plot(CLabeProbe_Neurons_V, CLabeProbe_Neurons_S, 'or');
    xlabel('Retinal Locaton');
    ylabel('Saccade Location');
    
    saveas(f,[netDir filesep 'CLayerProbe-summary.png']);
    
    close(f);
    
    analysisSummary = 0;
    
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

        saveas(f,[netDir filesep name '-summary-1.png']);
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

end