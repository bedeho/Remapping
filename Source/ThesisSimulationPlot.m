
%
%  ThesisSimulationPlot.m
%  Remapping
%
%  Created by Bedeho Mender on 25/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [stmCtrlFigure, remScatFig, remTraceScatFig] = ThesisSimulationPlot(experimentFolder)

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    if(nargin<1)
        experimentFolder = [EXPERIMENTS_FOLDER 'baseline-multiRF3-tuneC-to-Rpsi/C_to_R_psi=6/BlankNetwork/'];
    end
    
    color = [67,82,163]/255; % {[67,82,163]/255; [238,48,44]/255};
    
    %% Load stimuli control file
    stmCtrlFile = [experimentFolder filesep 'analysis-basic-StimuliControl.mat'];
    if(exist(stmCtrlFile, 'file')),
        
        stmCtrl = load(stmCtrlFile);
        StimuliControl_Result = stmCtrl.StimuliControl_Result;
        
        stmCtrlFigure = figure('Units','pixels','position', [1000 1000 620 300]);
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
            hBar = bar(x,hist(latency,x),1.0,'stacked','LineStyle','none');
            set(hBar,'FaceColor', color);
        end
        
        %xlim([0 maxLatency]);
        pbaspect([1 0.3 1]);
        box off;

        hXLabel = xlabel('Time (s)');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 14);
    
    end
    
    %% Load Remapping file
    remFile = [experimentFolder filesep 'analysis-basic-DuhamelRemapping.mat'];
    if(exist(remFile, 'file')),
        
        rem = load(remFile);
        DuhamelRemapping_Result = rem.DuhamelRemapping_Result;
        [remLatFig, remScatFig, indexFig] = remappingPlots({DuhamelRemapping_Result}, {color});

    end
    
    % Load Duhamel Trace file
    duhamelRemappingTraceFile = [experimentFolder filesep 'analysis-basic-DuhamelRemappingTrace.mat'];
    if(exist(duhamelRemappingTraceFile, 'file')),
        
        duhamelRemapping = load(duhamelRemappingTraceFile);
        DuhamelRemappingTrace_Result = duhamelRemapping.DuhamelRemappingTrace_Result;
        
        [remTraceLatFig, remTraceScatFig, indexTraceFig] = remappingPlots({DuhamelRemappingTrace_Result}, {color});
    end
    
    % Load Kusonoki File
    KusonokiFile = [experimentFolder filesep 'analysis-basic-Kusonoki.mat'];
    if(exist(KusonokiFile, 'file')),
        
        Kusonoki = load(KusonokiFile);
        kusonokiSTIMAlignedAnalysis = Kusonoki.kusonokiSTIMAlignedAnalysis;
        kusonokiSACCAlignedAnalysis = Kusonoki.kusonokiSACCAlignedAnalysis;
        
        %{
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
        %}
    end
    
    % Load CLayer Probe File
    CLayerProbeFile = [experimentFolder filesep 'analysis-basic-CLayerProbe.mat'];
    if(exist(CLayerProbeFile, 'file')),
        
        CLayerProbe = load(CLayerProbeFile);
        CLabeProbe_Neurons_S = CLayerProbe.CLabeProbe_Neurons_S;
        CLabeProbe_Neurons_V = CLayerProbe.CLabeProbe_Neurons_V;
        
        %{
        f = figure;
        plot(CLabeProbe_Neurons_V, CLabeProbe_Neurons_S, 'or');
        xlabel('Retinal Locaton');
        ylabel('Saccade Location');

        saveas(f,[netDir filesep 'CLayerProbe-summary.png']);

        close(f);
        %}
    end
    
    

end
