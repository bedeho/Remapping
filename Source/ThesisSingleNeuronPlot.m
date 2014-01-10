
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
        
    %% prewired - stim onset
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
    
    %% prewired - truncation
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
    
    %% prewired - remapping
    
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
    period      = R_BASE + (15);
    epoch       = 1;
    neuron      = R_BASE + (-20);
    stimuliName = 'STIM-basic-SaccadeControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/activity-basic-SaccadeControl.mat'];
    %}
    
    
    %% prewired - remapping
    
    % h_0 = -5& s = 15, --> r = -20.
    
    h = -5;
    s = 15;
    r = h - s;
    
    experiment  = 'baseline-onsettune';
    
    %{
    % Remapping
    period      = 1;
    epoch       = 1;
    neuron      = R_BASE + r;
    stimuliName = 'STIM-basic-DuhamelRemappingTrace';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/S_presaccadic_onset=0.07/TrainedNetwork/activity-basic-DuhamelRemappingTrace.mat'];
    %}
    
    %{
    % Stimuli control
    period      = R_BASE + h;
    epoch       = 1;
    neuron      = R_BASE + r;
    stimuliName = 'STIM-basic-StimuliControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/S_presaccadic_onset=0.07/TrainedNetwork/activity-basic-StimuliControl.mat'];
    %}
    
    
    % Saccade control
    period      = R_BASE + s;
    epoch       = 1;
    neuron      = R_BASE + h;
    stimuliName = 'STIM-basic-SaccadeControl';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/S_presaccadic_onset=0.07/TrainedNetwork/activity-basic-SaccadeControl.mat'];
    
    
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
    
end