
%
%  ThesisPopulationResponsePlot.m
%  Remapping
%
%  Created by Bedeho Mender on 24/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function ThesisPopulationResponsePlot()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    %{
    % Basic: Blank
    period      = 1;
    epoch       = 1;
    experiment  = 'test';
    stimuliFile = [EXPERIMENTS_FOLDER experiment filesep 'STIM-basic-DuhamelRemapping' filesep 'stim.mat'];
    activityFile = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork/activity-basic-DuhamelRemapping.mat'];
    %}
    
    %{
    % Basic: Prewired
    period      = 1;
    epoch       = 1;
    experiment  = 'prewired';
    stimuliFile = [EXPERIMENTS_FOLDER experiment filesep 'STIM-basic-DuhamelRemapping' filesep 'stim.mat'];
    activityFile = [EXPERIMENTS_FOLDER experiment filesep 'baseline/PrewiredNetwork/activity-basic-DuhamelRemapping.mat'];
    %}
        
    % Basic
    period      = 1;
    epoch       = 1;
    experiment  = 'baseline-onsettune';
    stimuliFile = [EXPERIMENTS_FOLDER experiment filesep 'STIM-basic-DuhamelRemappingTrace' filesep 'stim.mat'];
    activityFile = [EXPERIMENTS_FOLDER experiment filesep 'S_presaccadic_onset=0.07/TrainedNetwork/activity-basic-DuhamelRemappingTrace.mat'];
    
    
    %% ====================================================================================================================
    
    % Load input files
    disp('Loading input files...');
    
    activity = load(activityFile);
    V_firingrate = activity.V_firing_history(:, :, period, epoch);
    R_firingrate = activity.R_firing_history(:, :, period, epoch);
    S_firingrate = activity.S_firing_history(:, :, period, epoch);
    R_N = activity.R_N;
    S_N = activity.S_N;
    
    total = [R_firingrate; S_firingrate; V_firingrate;];
    
    stimuli  = load(stimuliFile);
    dt = stimuli.dt;
    R_eccentricity = stimuli.R_eccentricity;
    S_eccentricity = stimuli.S_eccentricity;
    
    s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);
    numTimeSteps = length(stimuli.stimuli{period}.eyePositionTrace);
    eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
    retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';
    
    % Setup axis ticks
    ticks = 1:(0.100/dt):numTimeSteps;
    for l=1:length(ticks),
        tickLabels{l} = num2str(stepToTime(ticks(l), dt));
    end
    
    % Plot
    h = figure('Units','pixels','position', [1000 1000 620 570]);
    
    % Response plot
    responseAxis = subplot(2,1,1);
    imagesc(total);
    colorbar
    
    % Stimuli plot
    stimAxis = subplot(2,1,2);
    hold on;
    plot(1:numTimeSteps, eyePositionTrace, 'r');
    plot(1:numTimeSteps, retinalTargetTraces, 'b');
    
    % We must to the subplots first, and posiiton later, some matlab quirck
    set(responseAxis, 'Units', 'Pixels', 'pos', [60 106 500 450]); % [left bottom width height]
    set(stimAxis, 'Units', 'Pixels', 'pos', [60 45 500 60]); % [left bottom width height]
    
    axes(responseAxis);
    hold on;
    set(gca, 'YDir', 'reverse');
    
    % Saccade bar
    total_N = 2*R_N+S_N;
    if ~isempty(s), plot([s s],[ones(total_N,1) total_N*ones(total_N,1)],'r'); end
    
    % Vertical bars
    plot([1 numTimeSteps], R_N + [1 1], 'w--');
    plot([1 numTimeSteps], S_N + R_N + [1 1], 'w--');
    
    % Pretty plot
    set(gca,'XTick', []);
    %axis tight;
    
    % Stimuli trace
    axes(stimAxis);
    hold on;

    hXLabel = xlabel('Time (s)');
    ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
    
    yTick = -30:15:R_eccentricity;
    xlim([1 numTimeSteps]);
    set(gca, 'YDir','reverse', 'XTick', ticks, 'XTickLabel', tickLabels, 'YTick', yTick);
    
    set(hXLabel, 'FontSize', 12);
    
    set(gca, 'FontSize', 10);
    set(gca, 'TickDir', 'out');
    
    % Fix ticks
    R_Deg = (-R_eccentricity:10:R_eccentricity);
    R_Ticks = R_Deg + R_eccentricity + 1;
    
    S_Deg = ((-S_eccentricity + 10):10:S_eccentricity);
    S_Ticks = S_Deg + S_eccentricity + 1 + (2*R_eccentricity + 1);
    
    V_Deg = ((-R_eccentricity + 10):10:R_eccentricity);
    V_Ticks = V_Deg + R_eccentricity + 1 + (2*S_eccentricity + 1) + (2*R_eccentricity + 1);
    
    Ticks = [R_Ticks S_Ticks V_Ticks];
    
    TickLabels = cell(1, length(Ticks));
    for i=1:length(Ticks),
        
        if(i <= length(R_Deg)),
            deg = R_Deg(i);
        elseif(i <= (length(R_Deg) + length(S_Deg))),
            deg = S_Deg(i - length(R_Deg));    
        else
            deg = V_Deg(i - length(S_Deg) - length(R_Deg));
        end
        
        TickLabels{i} = num2str(deg);
    end
    
    axes(responseAxis);
    set(gca,'TickDir', 'out','YTick', Ticks, 'YTickLabel', TickLabels);
    
    text(-30,78,'Remapping Neurons (deg)','Rotation',90);
    text(-30,150,'Saccade Neurons (deg)','Rotation',90);
    text(-30,225,'Visual Neurons (deg)','Rotation',90);

end
