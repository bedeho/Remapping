
%
%  ThesisSimulationPlot.m
%  Remapping
%
%  Created by Bedeho Mender on 25/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [stmCtrlFigure, remScatFig, remTraceScatFig, kusonokiSACCFigure, kusonokiSTIMFigure, CLayerProbeFigure] = ThesisSimulationPlot(experimentFolder)

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    if(nargin<1)
        experimentFolder = [EXPERIMENTS_FOLDER 'prewired/baseline/BlankNetwork/'];
    end
    
    color = [67,82,163]/255; % {[67,82,163]/255; [238,48,44]/255};
    
    %% Load stimuli control file
    stmCtrlFile = [experimentFolder filesep 'analysis-basic-StimuliControl.mat'];
    if(exist(stmCtrlFile, 'file')),
        
        stmCtrl = load(stmCtrlFile);
        StimuliControl_Result = stmCtrl.StimuliControl_Result;
        
        stmCtrlFigure = figure('Units','pixels','position', [1000 1000 500 300]);
        latency = 1000*[StimuliControl_Result(:).latency];
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
            %set(hBar,'FaceColor', color);
        end
        
        %xlim([0 maxLatency]);
        pbaspect([1 0.6 1]);
        axis tight

        hXLabel = xlabel('Time (ms)');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 16);
        set(gca, 'FontSize', 14);
    
    end
    
    %% Load Remapping file
    remFile = [experimentFolder filesep 'analysis-basic-DuhamelRemapping.mat'];
    if(exist(remFile, 'file')),
        
        rem = load(remFile);
        DuhamelRemapping_Result = rem.DuhamelRemapping_Result;
        [remLatFig, remScatFig, indexFig] = remappingPlots({DuhamelRemapping_Result}, {color});
    else
        remScatFig = 0; % set to garbage
    end
    
    % Load Duhamel Trace file
    duhamelRemappingTraceFile = [experimentFolder filesep 'analysis-basic-DuhamelRemappingTrace.mat'];
    if(exist(duhamelRemappingTraceFile, 'file')),
        
        duhamelRemapping = load(duhamelRemappingTraceFile);
        DuhamelRemappingTrace_Result = duhamelRemapping.DuhamelRemappingTrace_Result;
        
        [remTraceLatFig, remTraceScatFig, indexTraceFig] = remappingPlots({DuhamelRemappingTrace_Result}, {color});
    else
        remTraceScatFig = 0; % set to garbage
    end
    
    % Load Kusonoki File
    KusonokiFile = [experimentFolder filesep 'analysis-basic-Kusonoki.mat'];
    if(exist(KusonokiFile, 'file')),
        
        % Load analysis data
        Kusonoki = load(KusonokiFile);
        kusonokiSTIMAlignedAnalysis = Kusonoki.kusonokiSTIMAlignedAnalysis;
        kusonokiSACCAlignedAnalysis = Kusonoki.kusonokiSACCAlignedAnalysis;
        ticks = Kusonoki.ticks;
        
        % STIM aligned
        kusonokiSTIMFigure = figure('Units','pixels','position', [1000 1000 620 300]);

        arr_stim = kusonokiSTIMAlignedAnalysis;
        hold on;
        errorbar(ticks, [arr_stim(:).current_mean], [arr_stim(:).current_std],'-or');
        errorbar(ticks, [arr_stim(:).future_mean], [arr_stim(:).future_std],'-ob');
        hLenged = legend('Stimulus in Current RF','Stimulus in Future RF');
        legend boxoff;
        hXLabel = xlabel('Time from saccade onset to stimulus off (s)');
        hYLabel = ylabel({'Average response 0-350ms';'after stimulus onset'});
        set([hYLabel hXLabel hLenged], 'FontSize', 14);
        ylim([-0.1 1.1]);
        
        % SACC aligned
        kusonokiSACCFigure = figure('Units','pixels','position', [1000 1000 620 300]);
        arr_sacc = kusonokiSACCAlignedAnalysis;
        hold on;
        errorbar(ticks, [arr_sacc(:).current_mean], [arr_sacc(:).current_std],'-or');
        errorbar(ticks, [arr_sacc(:).future_mean], [arr_sacc(:).future_std],'-ob');
        hLenged = legend('Stimulus in Current RF','Stimulus in Future RF');
        legend boxoff;
        hXLabel = xlabel('Time from saccade onset to stimulus off (s)');
        hYLabel = ylabel({'Average response 0-300ms';'after saccade onset'});
        set([hYLabel hXLabel hLenged], 'FontSize', 14);
        ylim([-0.1 1.1]);
    else
        kusonokiSACCFigure = 0; % set to garbage
        kusonokiSTIMFigure = 0;
    end
    
    % Load CLayer Probe File
    CLayerProbeFile = [experimentFolder filesep 'analysis-basic-CLayerProbe.mat'];
    if(exist(CLayerProbeFile, 'file')),
        
        CLayerProbe = load(CLayerProbeFile);
        CLabeProbe_Neurons_S = CLayerProbe.CLabeProbe_Neurons_S;
        CLabeProbe_Neurons_V = CLayerProbe.CLabeProbe_Neurons_V;
        R_max = CLayerProbe.R_max;
        S_max = CLayerProbe.S_max;
        
        CLayerProbeFigure = figure('Units','pixels','position', [1000 1000 420 300]);
        plot(CLabeProbe_Neurons_V, CLabeProbe_Neurons_S, 'or');
        hXLabel = xlabel('Retinal Locaton (deg)');
        hYLabel = ylabel('Saccade (deg)');
        set([hYLabel hXLabel], 'FontSize', 14);
        
        xlim(1.1*[-R_max R_max]);
        ylim(1.1*[-S_max S_max]);
        pbaspect([2*R_max 2*S_max 1]);
    else
        CLayerProbeFigure = 0; % set to garbage
    end
    
    

end
