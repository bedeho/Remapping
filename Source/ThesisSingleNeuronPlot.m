
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
    %global STIMULI_FOLDER;
    global EXPERIMENTS_FOLDER;
    
    % activityFiles{1} = [EXPERIMENTS_FOLDER 'test/baseline/BlankNetwork/activity-basic-DuhamelRemapping.mat'];
    
    % Basic
    period      = 1;
    epoch       = 1;
    neuron      = 46;
    experiment  = 'test';
    stimuliName = 'STIM-basic-DuhamelRemapping';
    activityFiles{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork/activity-basic-DuhamelRemapping.mat'];
    activityFiles{2} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/activity-basic-DuhamelRemapping.mat'];
    colors      = {'-b'; '--k'};
    legends     = {'Untrained'; 'Trained'};
    
    %%
    
    % Load input files
    disp('Loading input files...');
    
    %stimuli  = load([STIMULI_FOLDER stimuliName filesep 'stim.mat']);
    stimuli  = load([EXPERIMENTS_FOLDER experiment filesep stimuliName filesep 'stim.mat']);
    dt = stimuli.dt;

    s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);
    numTimeSteps = length(stimuli.stimuli{period}.eyePositionTrace);

    % Plot
    h = figure('Units','pixels','position', [1000 1000 600 250]);
    
    % Setup axis ticks
    ticks = 1:(0.100/dt):numTimeSteps;
    for l=1:length(ticks),
        tickLabels{l} = num2str(stepToTime(ticks(l), dt));
    end
    
    % Iterate activity files and plot
    responseAxis = subplot(2,1,1);
    set(responseAxis, 'Units', 'Pixels', 'pos', [60 120 500 110]); % [left bottom width height]
    hold on;
    
    for i=1:length(activityFiles),
        
        activity = load(activityFiles{i});
        
        % Load
        responseTrace = activity.R_firing_history(neuron, :, period, epoch);
        
        % Plot
        plot(1:numTimeSteps, responseTrace, colors{i});
        
    end
    
    % Pretty up
    hXLabel = xlabel('Time (s)');
    hYLabel = ylabel('Firing Rate');
    ylim([-0.05 max(1, max(responseTrace))]);
    xlim([1 numTimeSteps]);
    legend(legends);
    legend('boxoff');

    % Saccade times 
    if ~isempty(s), plot([s s],[-0.05 1],'r'); end 
    
    set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
    set(gca,'YTick', [0 1]);

    set([hYLabel hXLabel], 'FontSize', 14);
    
    set([gca], 'FontSize', 14);
    set( gca, 'TickDir', 'out' );
    

    % Add bottom traces
    eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
    retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';

    traceAxis = subplot(2,1,2);
    set(traceAxis, 'Units', 'Pixels', 'pos', [60 35 500 40]); % [left bottom width height]
    hold on;
    
    plot(1:numTimeSteps, eyePositionTrace, 'r');

    if(~isempty(retinalTargetTraces)),
        plot(0:(numTimeSteps-1), retinalTargetTraces, 'b');       
    end

    
    ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
    xlim([1 numTimeSteps]);
    set(gca,'YDir','reverse');
    set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
    axis off
    
end