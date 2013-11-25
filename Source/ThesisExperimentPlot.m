
%
%  ThesisExperimentPlot.m
%  Remapping
%
%  Created by Bedeho Mender on 25/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function ThesisExperimentPlot(experimentFolder)

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    if(nargin<1)
        experimentFolder = [EXPERIMENTS_FOLDER experiment 'test/baseline/BlankNetwork/'];
    end
    
    % Load stimuli control file
    stmCtrlFile = [experimentFolder filesep 'analysis-StimuliControl.mat'];
    if(exist(stmCtrlFile)),
        stmCtrl = load(stmCtrlFile);
        StimuliControl_Result = stmCtrl.StimuliControl_Result;
    end
    
    % Load saccade control file
    sacCtrlFile = [experimentFolder filesep 'analysis-SaccadeControl.mat'];
    if(exist(sacCtrlFile)),
        sacCtrl = load(sacCtrlFile);
        SaccadeControl_Result = sacCtrl.SaccadeControl_Result;
    end
    
    % Load Remapping file
    remFile = [experimentFolder filesep 'analysis-DuhamelRemapping.mat'];
    if(exist(remFile)),
        rem = load(remFile);
        DuhamelRemapping_Result = rem.DuhamelRemapping_Result;
    end
    
    % Load Duhamel Trace file
    [experimentFolder filesep 'analysis-DuhamelRemappingTrace.mat']
    
    
    DuhamelRemappingTrace_Result = save( , 'DuhamelRemappingTrace_Result');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    Kusonoki_Result = save([experimentFolder filesep 'analysis-Kusonoki.mat'] )
    
    , 'kusonokiSTIMAlignedAnalysis', 'kusonokiSACCAlignedAnalysis'
    
    
    [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activity, stimuli);
            
            ;
    
    
                        [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity,  stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'CLabeProbe_Neurons_S', 'CLabeProbe_Neurons_V');
            

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

end
