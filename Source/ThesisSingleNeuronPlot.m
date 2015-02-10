
%
%  ThesisSingleNeuronPlot.m
%  SMI
%
%  Created by Bedeho Mender on 23/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function ThesisSingleNeuronPlot()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    R_BASE = 46;
        
    %% test
    %{
    experiment = 'test';
    
    period      = R_BASE + (-20);
    epoch       = 5;
    neuron      = period;
    %stimuliName = 'STIM-basic-StimuliControl';
    %activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-StimuliControl.mat'];
    stimuliName = 'STIM-basic-Kusonoki';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-Kusonoki.mat'];

    colors{1}   = [0,0,255]/255; % [67,82,163]/255;
    legends{1}  = '';
    %}
    
    %{
    experiment  = 'test';
    
    % Remapping
    period      = 5;
    epoch       = 1;
    neuron      = 26;
    stimuliName = 'STIM-basic-Kusonoki';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-Kusonoki.mat'];
    %}
    
    %% stim onset (prewired)
    %{
    experiment = 'prewired';
    
    period      = R_BASE + (-20);
    epoch       = 1;
    neuron      = period;
    stimuliName = 'STIM-basic-StimuliControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-StimuliControl.mat'];

    colors{1}   = [0,0,255]/255; % [67,82,163]/255;
    legends{1}  = '';
    %}
    
    %% truncation (prewired)
    %{
    experiment = 'prewired';
    
    period      = 1;
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-DuhamelTruncation';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-DuhamelTruncation.mat'];

    colors{1}   = [0,0,255]/255; % [67,82,163]/255;
    legends{1}  = '';
    %}
    
    %% remapping (prewired)
    %{
   
    % h_0 = -5& s = 15, --> r = -20.
    experiment  = 'prewired';
    
    % Remapping
    %{
    period      = 1;
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-DuhamelRemappingTrace';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-DuhamelRemappingTrace.mat'];
    %}
    
    % Remapping 2
    %{
    period      = 1;
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-DuhamelRemapping';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-DuhamelRemapping.mat'];
    %}
    
    % Stimuli control
    %{
    period      = R_BASE + (-5);
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-StimuliControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-StimuliControl.mat'];
    %}
    
    % Saccade control
    %{
    period      = 46; % not contigous due to minimal saccade restriction: R_BASE + (15);
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-SaccadeControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-SaccadeControl.mat'];
    %}
    
    %}
    
    %% remapping
    
    % h_0 = -5& s = 15, --> r = -20.
    
    h = -5;
    s = 15;
    r = h - s;
    
    experiment  = 'classic'; % sprattling_visual_learning_bigepoch20-connectivitfix
    
    %{ 
    % Remapping
    period      = 1;
    epoch       = 1;
    neuron      = R_BASE + r;
    stimuliName = 'STIM-basic-DuhamelRemappingTrace';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-DuhamelRemappingTrace.mat'];
    %}
    
    %{
    % Stimuli control
    period      = R_BASE + h;
    epoch       = 1;
    neuron      = R_BASE + r;
    stimuliName = 'STIM-basic-StimuliControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-StimuliControl.mat'];
    %}
    
    % Saccade control
    
    period      = 46;
    epoch       = 1;
    neuron      = R_BASE + h;
    stimuliName = 'STIM-basic-SaccadeControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-SaccadeControl.mat'];
    
    
    
    
    %% KUSONKI
    
    %{
    
    % h_0 = -5& s = 15, --> r = -20.
    
    h = -5;
    s = 15;
    r = h - s;
    
    %experiment  = 'baseline-onsettune';
    experiment  = 'classic';
    
    % CRF periods: [1,6,8,13]
    % FRF periods: 13 + [1,6,8,13]
    
    % Remapping
    period      = 13 + 13;
    epoch       = 1;
    neuron      = R_BASE + r;
    stimuliName = 'STIM-basic-Kusonoki';
    %activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/S_presaccadic_onset=0.07/TrainedNetwork/activity-basic-Kusonoki.mat'];
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-Kusonoki.mat'];
    
    %}
    colors{1}   = [0,0,255]/255; % [67,82,163]/255;
    legends{1}  = '';
    
    
    
    % =======================================
    
    % Load input files
    disp('Loading input files...');
    
    stimuli  = load([EXPERIMENTS_FOLDER experiment filesep stimuliName filesep 'stim.mat']);
    dt = stimuli.dt;

    s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);
    numTimeSteps = length(stimuli.stimuli{period}.eyePositionTrace);
    
    width = 600; % classic = 600
    panel_width = width - 100;

    % Plot
    h = figure('Units','pixels','position', [1000 1000 width 250]); % [1000 1000 600 250]
    
    % Setup axis ticks
    ticks = 1:(0.100/dt):numTimeSteps;
    for l=1:length(ticks),
        tickLabels{l} = num2str(stepToTime(ticks(l), dt));
    end
    
    % Iterate activity files and plot
    responseAxis = subplot(2,1,1);
    set(responseAxis, 'Units', 'Pixels', 'pos', [60 120 panel_width 110]); % [left bottom width height] = classic [60 120 500 110]
    hold on;
    
    % plot rectangle
    stimBegin = stimuli.stimuli{period}.stimOnsetTimes;
    
    
    if(~isempty(stimBegin)),

        recordingWindow = 0.300;
        
        stimColor = 1.2*[0, 164, 103]/255;
        stimWindowBegin = stimBegin + 0.050;
        pos_box = [timeToTimeStep(stimWindowBegin, dt), 0.001, timeToTimeStep(recordingWindow, dt),  1];
        
        %rectangle('Position', pos_box,'FaceColor', stimColor,'EdgeColor', stimColor);
        
        %[x1, x2, x3, x4] = pos_box_to_X(pos_box);
        %[y1, y2, y3, y4] = pos_box_to_Y(pos_box);
        
        %p = patch([x1, x2, x3, x4], [y1, y2, y3, y4], stimColor);
        %set(p, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', stimColor, 'EdgeColor', stimColor);
        
        %{
        saccColor = [242, 240, 94]/255;
        saccBegin = stimuli.stimuli{period}.saccadeTimes;
        pos_box = [timeToTimeStep(saccBegin, dt), 0.001, timeToTimeStep(recordingWindow, dt),  1];
        rectangle('Position', pos_box,'FaceColor', saccColor,'EdgeColor', saccColor);
        %}
        
        %X = pos_box_to_X(pos_box);
        %Y = pos_box_to_Y(pos_box);
        
        %p = patch(X, Y, stimColor);
        %set(p,'FaceAlpha',0.5);
    end
    
    
    for i=1:length(activityFiles),
        
        activity = load(activityFiles{i});
        
        % Load
        responseTrace = activity.R_firing_history(neuron, :, period, epoch);
        
        % Plot
        plot(1:numTimeSteps, responseTrace, 'Color', colors{i});
        
    end
    
    % Pretty up
    hXLabel = xlabel('Time (s)');
    hYLabel = ylabel('Firing Rate');
    ylim([-0.05 max(1, max(responseTrace))]);
    xlim([1 numTimeSteps]);
    
    if(length(legends) > 1)
        legend(legends);
        legend boxoff;
    end

    % Saccade times 
    if ~isempty(s), plot([s s],[-0.05 1],'r'); end 
    
    set(gca,'XTick', ticks, 'XTickLabel', tickLabels, 'YTick', [0 1], 'TickDir', 'out', 'FontSize', 14);
    set([hYLabel hXLabel], 'FontSize', 14);

    % Add bottom traces
    eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
    retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';

    traceAxis = subplot(2,1,2);
    set(traceAxis, 'Units', 'Pixels', 'pos', [60 35 panel_width 40]); % [left bottom width height] = [60 35 500 40]
    hold on;
    
    plot(1:numTimeSteps, eyePositionTrace, 'r');

    if(~isempty(retinalTargetTraces)),
        plot(1:numTimeSteps, retinalTargetTraces, 'b');       
    end
    
    ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
    xlim([1 numTimeSteps]);
    set(gca,'YDir','reverse');
    set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
    axis off
    
    function [x1, x2, x3, x4] = pos_box_to_X(box),
        
        x1 = box(1);
        x2 = x1;
        x3 = x1 + box(3);
        x4 = x3;
       
    end
    
    function [y1, y2, y3, y4] = pos_box_to_Y(box),
       
        y1 = box(2);
        y2 = y1 + box(4);
        y3 = y2;
        y4 = y1;
        
    end
end