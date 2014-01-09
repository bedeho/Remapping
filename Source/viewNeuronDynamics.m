
%
%  viewNeuronDynamics.m
%  Remapping
%
%  Created by Bedeho Mender on 26/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function viewNeuronDynamics(activityFile, stimuliFile, networkFile, CLayerProbleFile)

    % Import global variables
    declareGlobalVars();
    
    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-KusonokiTesting.mat';
    end

    % CLayerProbe
    %[pathstr, name, ext] = fileparts(activityFile);
    %CLayerProbleFile = [pathstr filesep 'analysis-basic-CLayerProbe.mat']
    
    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load([stimuliFile filesep 'stim.mat']);
    network = load(networkFile);
    stimuliType = stimuli.stimulitype;
    
    if nargin == 4 && exist(CLayerProbleFile), %% HACK
        CLayerProbleFileAnalysis = load(CLayerProbleFile);
        
        S = CLayerProbleFileAnalysis.CLabeProbe_Neurons_S;
        V = CLayerProbleFileAnalysis.CLabeProbe_Neurons_V;
        R_max = CLayerProbleFileAnalysis.R_max;
        S_max = CLayerProbleFileAnalysis.S_max;
    end
    
    % Set parameters
    R_N = activity.R_N;
    C_N = activity.C_N;
    dt = stimuli.dt;
    
    % Make figure
    figure('name',stimuliType,'Position', [100, 100, 1049, 895]);
    
    % Adding controls
    if activity.numEpochs > 1
        
        menu = '1';
        
        for e=2:activity.numEpochs,
            menu = [menu '|' num2str(e)];
        end
        
        %uicontrol('Style', 'slider', 'Min', 1, 'Max', activity.numEpochs, 'Value', 1, 'Position', [20 850 120 20], 'Callback', {@epoch_callback});
        
        uicontrol('Style', 'text', 'String', 'Epoch', 'Position', [20 865 120 20]);
        uicontrol('Style', 'popup', 'String', menu, 'Position', [20 850 120 20], 'Callback', @epoch_callback);
    end
    
    if activity.numPeriods > 1
        %uicontrol('Style', 'slider', 'Min', 1, 'Max', activity.numPeriods, 'Value', 1, 'Position', [200 850 120 20], 'Callback', {@period_callback});
        
        menu = '1';
        
        for p=2:activity.numPeriods,
            menu = [menu '|' num2str(p)];
        end
        
        uicontrol('Style', 'text', 'String', 'Period', 'Position', [20 355 100 50]);
        uicontrol('Style', 'popup', 'String', menu, 'Position', [20 340 100 50], 'Callback', @period_callback);
    end
    
    
    function epoch_callback(hObj,event,ax)
        
        epoch = get(hObj,'Value');
        display();
    end 
    
    function period_callback(hObj,event,ax)
        
        period = get(hObj,'Value');
        display();
    end

    % Setup global vars
    period = 1
    epoch = 1
    imgV = [];
    imgR = [];
    imgS = [];
    imgC = [];
    
    % Fix R axes ticks
    if(isfield(stimuli,'R_eccentricity'))
        R_eccentricity = stimuli.R_eccentricity;
    else
        R_eccentricity = 45;
    end
    R_preferences = -R_eccentricity:1:R_eccentricity;
    
    rTicks = 1:length(R_preferences);
    rdist = 15;
    rTicks = rTicks(1:rdist:end);
    rLabels = R_preferences(1:rdist:end);
    rCellLabels = cell(1,length(rLabels));
    for i=1:length(rLabels),
      rCellLabels{i} = num2str(rLabels(i));
    end
    
    % Fix S axes ticks
    if(isfield(stimuli,'S_eccentricity'))
        S_eccentricity = stimuli.S_eccentricity;
    else
        S_eccentricity = 30;
    end
    
    S_preferences = -S_eccentricity:1:S_eccentricity;
    
    sTicks = 1:length(S_preferences);
    sdist = 5;
    sTicks = sTicks(1:sdist:end);
    sLabels = S_preferences(1:sdist:end);
    sCellLabels = cell(1,length(sLabels));
    for i=1:length(sLabels),
      sCellLabels{i} = num2str(sLabels(i));
    end
    
    % Do first plot
    display();

    function display()
        
        %{
        % Get desired epoch or quit
        if activity.numEpochs > 1,
            str = input(['Epoch # (1, ' num2str(activity.numEpochs) '), or 0 to quit: '],'s');
            epoch = str2num(str);
            
            % Quit if we are done
            if epoch == 0,
                return;
            end
        else
            epoch = 1;
        end
        
        % Get desired period
        str = input(['Period # (1, ' num2str(activity.numPeriods) '): '],'s');
        period = str2num(str);
        %}
        
        stimuli.stimuli{period}
        
        % Load
        V_firingrate = activity.V_firing_history(:, :, period, epoch);
        R_firingrate = activity.R_firing_history(:, :, period, epoch);
        S_firingrate = activity.S_firing_history(:, :, period, epoch);
        V_activation = activity.V_activation_history(:, :, period, epoch);
        R_activation = activity.R_activation_history(:, :, period, epoch);
        S_activation = activity.S_activation_history(:, :, period, epoch);
        
        extra = activity.extra_history(:, :, period, epoch);
        
        % If C is empty, just fill with blank
        if ~isempty(activity.C_firing_history),
            C_firingrate = activity.C_firing_history(:, :, period, epoch);
            C_activation = activity.C_activation_history(:, :, period, epoch);
            
            max(max(C_activation))
            
            includeC_layer = 1;
            %numRows = 2+5;
            
        else
            C_firingrate = zeros(C_N, activity.numPeriods);
            C_activation = zeros(C_N, activity.numPeriods);
            
            includeC_layer = 0;
            %numRows = 5;
        end
        
        numRows = 2+5;
        
        % Plot
        s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);
        numTimeSteps = length(stimuli.stimuli{period}.eyePositionTrace);
        
        ticks = 1:(0.100/dt):numTimeSteps;
        %ticks(2:end) = ticks(2:end) + 1;
        for l=1:length(ticks),
            tickLabels{l} = num2str(stepToTime(ticks(l), dt)); % (l-1)*0.100
        end
        
        subplot(numRows,2,1);
        imgR = imagesc(extra);
        hold on;colorbar;
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('Extra');
        set(imgR, 'ButtonDownFcn', {@singleUnitCallBack, 'E'}); % Setup callback
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);

        %{
        subplot(numRows,2,2);
        imagesc(E_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('E Firing');
        %}
        
        subplot(numRows,2,3);
        imgV = imagesc(V_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('V Firing');
        set(imgV, 'ButtonDownFcn', {@singleUnitCallBack, 'V'}); % Setup callback
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        subplot(numRows,2,4);
        imagesc(V_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('V Activation');
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        subplot(numRows,2,5);
        imgR = imagesc(R_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Firing');
        set(imgR, 'ButtonDownFcn', {@singleUnitCallBack, 'R'}); % Setup callback
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        subplot(numRows,2,6);
        imagesc(R_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Activation');
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        subplot(numRows,2,7);
        imgS = imagesc(S_firingrate);
        colorbar
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('S Firing');
        set(imgS, 'ButtonDownFcn', {@singleUnitCallBack, 'S'}); % Setup callback
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        subplot(numRows,2,8);
        imagesc(S_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('S Activation');
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        
        %if(includeC_layer),
            subplot(numRows,2,[9 11]);
            imgC = imagesc(C_firingrate);
            hold on;colorbar
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Firing');
            set(imgC, 'ButtonDownFcn', {@singleUnitCallBack, 'C'}); % Setup callback
            set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
            
            subplot(numRows,2,[10 12]);
            imagesc(C_activation);
            hold on;colorbar; caxis([min(min(C_activation)) max(max(C_activation))]);
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Actiation');
            set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
            
            nextplot = 13;
        %else
        %    nextplot = 9;
        %end
        
        % Add bottom traces
        eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
        retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';

        subplot(numRows,2,nextplot);
        cla
        plot(1:numTimeSteps, eyePositionTrace, 'r');
        hold on;
        
        if(~isempty(retinalTargetTraces)),
        
            plot(1:numTimeSteps, retinalTargetTraces, 'b');
            legend({'Eye Position','Stimuli Retinal Locations'});
        else
            
            legend({'Eye Position'});            
        end
        
        xlabel(['Time step (s)']);
        ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
        xlim([1 numTimeSteps]);
        set(gca,'YDir','reverse');
        set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
        %{
        subplot(numRows,2,nextplot+1);
        cla
        plot(0:(numTimeSteps-1), eyePositionTrace, 'r');
        hold on;
        plot(0:(numTimeSteps-1), retinalTargetTraces, 'b');
        xlabel(['Time step (dt =' num2str(dt) ')']);
        legend({'Eye Position','Stimuli Retinal Locations'});
        ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
        xlim([0 (numTimeSteps-1)]);
        set(gca,'YDir','reverse');
        %}

        function singleUnitCallBack(varargin)
            %
            % Extract region,row,col
            region = varargin{3}

            % Pick neuron
            x = inputdlg('Neuron #:', 'Sample', [1 50]);
            neuron = str2num(cell2mat(x))

            % single left  click => 'SelectionType' = 'normal'
            % single right click => 'SelectionType' = 'alt'
            % double right click => 'SelectionType' = 'open'
            clickType = get(gcf,'SelectionType');
            
            % Check which pool we are looking atr
            if(strcmp(region,'V')),
                responseTrace = activity.V_firing_history(neuron, :, period, epoch);
            elseif(strcmp(region,'E')),
                responseTrace = activity.extra_history(neuron, :, period, epoch);
            elseif(strcmp(region,'R')),
                responseTrace = activity.R_firing_history(neuron, :, period, epoch);
            elseif(strcmp(region,'S')),
                responseTrace = activity.S_firing_history(neuron, :, period, epoch);
            elseif(strcmp(region,'C')),
                
                if(includeC_layer),
                    responseTrace = activity.C_firing_history(neuron, :, period, epoch);
                else
                    responseTrace = ones(1,numTimeSteps);
                end
            end
            
            if(strcmp(clickType, 'normal')), % response trace

                %
                % response trace
                %
                
                figure('Units','normalized','position',[.1 .1 .22 .1]);
                hold on;

                eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
                retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';
            
                targetOffIntervals = stimuli.stimuli{period}.targetOffIntervals{1};
                [numOffPeriods,~] = size(targetOffIntervals);
                
                lastStimOnsetTime = [];
                for t=1:numOffPeriods,

                    if(isempty(lastStimOnsetTime))
                        lastStimOnsetTime = targetOffIntervals(t,2);
                    else
                        % plot rectangle
                        stimDuration = targetOffIntervals(t,1)-lastStimOnsetTime;
                        rectangle('Position', [timeToTimeStep(lastStimOnsetTime, dt), 0.001, timeToTimeStep(stimDuration, dt),  1],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);

                        lastStimOnsetTime = [];
                    end
                end

                if(~isempty(lastStimOnsetTime))
                    rectangle('Position', [timeToTimeStep(lastStimOnsetTime, dt), 0.001, timeToTimeStep(numTimeSteps-lastStimOnsetTime, dt), 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
                end

                plot(1:numTimeSteps, responseTrace, 'b');
                if ~isempty(s), plot([s s],[-0.05 1],'r'); end % Saccade times
                hXLabel = xlabel('Time (s)');
                hYLabel = ylabel('Firing Rate');
                ylim([-0.05 max(1, max(responseTrace))]);
                xlim([1 numTimeSteps]);
                
                set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
                set(gca,'YTick', [0 1]);

                set([hYLabel ], 'FontSize', 14);
                set([hXLabel], 'FontSize', 16);
                set([gca], 'FontSize', 14);
                set( gca, 'TickDir', 'out' );
                
                %
                % Stimuli trace
                %
                
                figure('Units','normalized','position',[.1 .1 .22 .1]);
                hold on;

                eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
                retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';
            
                plot(1:numTimeSteps, eyePositionTrace, 'r');

                hXLabel = xlabel('Time (s)');
                hYLabel = ylabel('Location (deg)');
                ylim([-45 45]);
                xlim([1 numTimeSteps]);
                
                if(~isempty(retinalTargetTraces)),
                    plot(1:numTimeSteps, retinalTargetTraces, 'b');
                    legend({'Eye','Stimuli'});
                else

                    legend({'Eye'});            
                end
                
                
                legend('Location','NorthWest')
                legend boxoff
                set(gca,'XTick', ticks, 'XTickLabel', tickLabels);
                set(hYLabel, 'FontSize', 14);
                set(hXLabel, 'FontSize', 14);
                set(gca, 'FontSize', 10, 'TickDir', 'out' ,'YDir','reverse');
                
            else % weight vectors!
                
                if(strcmp(region,'V')),
                    
                    figure('Units','pixels','position', [1000 1000 420 300]);
                    hold on;
                    
                    synapticEfferents = network.V_to_C_weights(:, neuron);
                    maxWeight = max(synapticEfferents);
                    num = length(synapticEfferents);

                    for i=1:num,
                        weight = synapticEfferents(i);
                        plot(V(i),S(i),'o','Color',[1, 1 -  weight/maxWeight, 1 -  weight/maxWeight]);
                    end
                    
                    hXLabel = xlabel('Retinal Locaton (deg)');
                    hYLabel = ylabel('Saccade (deg)');
                    set([hYLabel hXLabel], 'FontSize', 14);

                    box on;
                    
                    xlim(1.1*[-R_max R_max]);
                    ylim(1.1*[-S_max S_max]);
                    pbaspect([2*R_max 2*S_max 1]);
                    
                elseif(strcmp(region,'R')),
                    
                    figure('Units','pixels','position', [1000 1000 420 230]);
                    hold on;
                    
                    synapticAfferents = network.C_to_R_weights(neuron, :);
                    maxWeight = max(synapticAfferents);
                    num = length(synapticAfferents);
                    
                    for i=1:num,
                        weight = synapticAfferents(i);
                        plot(V(i),S(i),'o','Color',[1, 1 -  weight/maxWeight, 1 -  weight/maxWeight]);
                    end
                    
                    hXLabel = xlabel('Retinal Locaton (deg)');
                    hYLabel = ylabel('Saccade (deg)');
                    set([hYLabel hXLabel], 'FontSize', 14);

                    box on;
                    
                    xlim(1.1*[-R_max R_max]);
                    ylim(1.1*[-S_max S_max]);
                    pbaspect([2*R_max 2*S_max 1]);
                    
                    colorbar;
                    caxis([0 maxWeight]);
                    resolution = 100;
                    map = [ones(101,1) (1:-1/resolution:0)' (1:-1/resolution:0)'];
                    colormap(map);
                    
                    % Diag
                    disp('C afferents > 90% of max');
                    find(synapticAfferents > maxWeight*0.9)

                elseif(strcmp(region,'S')),
                    
                    figure('Units','pixels','position', [1000 1000 420 300]);
                    hold on;
                    
                    synapticEfferents = network.S_to_C_weights(:, neuron);
                    maxWeight = max(synapticEfferents);
                    num = length(synapticEfferents);

                    for i=1:num,
                        weight = synapticEfferents(i);
                        plot(V(i),S(i),'o','Color',[1, 1 -  weight/maxWeight, 1 -  weight/maxWeight]);
                    end
                    
                    hXLabel = xlabel('Retinal Locaton (deg)');
                    hYLabel = ylabel('Saccade (deg)');
                    set([hYLabel hXLabel], 'FontSize', 14);
                    set(gca, 'FontSize', 12);
                    
                    box on;
                    
                    xlim(1.1*[-R_max R_max]);
                    ylim(1.1*[-S_max S_max]);
                    pbaspect([2*R_max 2*S_max 1]);
                    
                elseif(strcmp(region,'C')),
                    
                    C_to_R_weightvector = network.C_to_R_weights(:, neuron);
                    S_to_C_afferents =  network.S_to_C_weights(neuron, :);
                    V_to_C_afferents = network.V_to_C_weights(neuron, :);

                    % 'S->C'
                    figure('Units','pixels','position', [1000 1000 420 240]);
                    plot(S_to_C_afferents);
                    axis tight;
                    ylim([0 max(S_to_C_afferents)]);
                    
                    hXLabel = xlabel('Saccade (deg)');
                    hYLabel = ylabel('Synaptic weight');
                    set([hYLabel hXLabel], 'FontSize', 14);
                    set(gca,'XTick', sTicks, 'XTickLabel', sCellLabels);
                    pbaspect([1, 30/length(S_preferences), 1]);
                    
                    % 'V->C'
                    figure('Units','pixels','position', [1000 1000 420 180]);
                    plot(V_to_C_afferents);
                    axis tight;
                    ylim([0 max(V_to_C_afferents)]);
                    
                    hXLabel = xlabel('Retinal Location (deg)');
                    hYLabel = ylabel('Synaptic weight');
                    set([hYLabel hXLabel], 'FontSize', 14);
                    set(gca,'XTick', rTicks, 'XTickLabel', rCellLabels);
                    pbaspect([1, 30/length(R_preferences), 1]);
                    
                    
                    % 'C->R'
                    figure('Units','pixels','position', [1000 1000 420 180]);
                    plot(C_to_R_weightvector);
                    axis tight;
                    ylim([0 max(C_to_R_weightvector)]);
                    
                    hXLabel = xlabel('Retinal Location (deg)');
                    hYLabel = ylabel('Synaptic weight');
                    set([hYLabel hXLabel], 'FontSize', 14);
                    set(gca,'XTick', rTicks, 'XTickLabel', rCellLabels);
                    pbaspect([1, 30/length(R_preferences), 1]);
                    
                end

            end
        end

    end


end

