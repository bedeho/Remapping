
%
%  viewNeuronDynamics.m
%  Remapping
%
%  Created by Bedeho Mender on 26/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function viewNeuronDynamics(activityFile, stimuliName, networkFile)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;
    
    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-KusonokiTesting.mat';
        stimuliName     = 'basic-KusonokiTesting';
    end

    % CLayerProbe
    [pathstr, name, ext] = fileparts(activityFile);
    CLayerProbleFile = [pathstr filesep 'analysis-basic-CLayerProbe.mat']
    
    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load([STIMULI_FOLDER stimuliName filesep 'stim.mat']);
    network = load(networkFile);
    CLayerProbleFileAnalysis = load(CLayerProbleFile);
    
    % Set parameters
    R_N = activity.R_N;
    C_N = activity.C_N;
    dt = stimuli.dt;
    
    % Make figure
    figure('name',stimuliName,'Position', [100, 100, 1049, 895]);
    
    % Adding controls
    if activity.numEpochs > 1
        uicontrol('Style', 'slider', 'Min', 1, 'Max', activity.numEpochs, 'Value', 1, 'Position', [20 850 120 20], 'Callback', {@epoch_callback});
    end
    
    if activity.numPeriods > 1
        %uicontrol('Style', 'slider', 'Min', 1, 'Max', activity.numPeriods, 'Value', 1, 'Position', [200 850 120 20], 'Callback', {@period_callback});
        
        menu = '1';
        
        for p=2:activity.numPeriods,
            menu = [menu '|' num2str(p)];
        end
        
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
        E_firingrate = activity.E_firing_history(:, :, period, epoch);
        V_firingrate = activity.V_firing_history(:, :, period, epoch);
        R_firingrate = activity.R_firing_history(:, :, period, epoch);
        S_firingrate = activity.S_firing_history(:, :, period, epoch);
        
        E_activation = activity.E_activation_history(:, :, period, epoch);
        V_activation = activity.V_activation_history(:, :, period, epoch);
        R_activation = activity.R_activation_history(:, :, period, epoch);
        S_activation = activity.S_activation_history(:, :, period, epoch);
        
        extra = activity.extra_history(:, :, period, epoch);
        
        % If C is empty, just fill with blank
        if ~isempty(activity.C_firing_history),
            C_firingrate = activity.C_firing_history(:, :, period, epoch);
            C_activation = activity.C_activation_history(:, :, period, epoch);
            
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

        subplot(numRows,2,1);
        imagesc(extra);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('Extra');

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

        
        subplot(numRows,2,4);
        imagesc(V_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('V Activation');
        
        
        subplot(numRows,2,5);
        imgR = imagesc(R_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Firing');
        set(imgR, 'ButtonDownFcn', {@singleUnitCallBack, 'R'}); % Setup callback

        subplot(numRows,2,6);
        imagesc(R_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Activation');

        subplot(numRows,2,7);
        imgS = imagesc(S_firingrate);
        colorbar
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('S Firing');
        set(imgS, 'ButtonDownFcn', {@singleUnitCallBack, 'S'}); % Setup callback

        subplot(numRows,2,8);
        imagesc(S_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('S Activation');

        %if(includeC_layer),
            subplot(numRows,2,[9 11]);
            imgC = imagesc(C_firingrate);
            hold on;colorbar
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Firing');
            set(imgC, 'ButtonDownFcn', {@singleUnitCallBack, 'C'}); % Setup callback

            subplot(numRows,2,[10 12]);
            imagesc(C_activation);
            hold on;colorbar
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Actiation');
            
            nextplot = 13;
        %else
        %    nextplot = 9;
        %end
        
        % Add bottom traces
        eyePositionTrace = stimuli.stimuli{period}.eyePositionTrace;
        retinalTargetTraces = stimuli.stimuli{period}.retinalTargetTraces';

        subplot(numRows,2,nextplot);
        cla
        plot(0:(numTimeSteps-1), eyePositionTrace, 'r');
        hold on;
        
        if(~isempty(retinalTargetTraces)),
        
            plot(0:(numTimeSteps-1), retinalTargetTraces, 'b');
            legend({'Eye Position','Stimuli Retinal Locations'});
        else
            
            legend({'Eye Position'});            
        end
        
        xlabel(['Time step (dt =' num2str(dt) ')']);
        ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
        xlim([0 (numTimeSteps-1)]);
        set(gca,'YDir','reverse');
        
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
            neuron = str2num(cell2mat(x));

            % single left  click => 'SelectionType' = 'normal'
            % single right click => 'SelectionType' = 'alt'
            % double right click => 'SelectionType' = 'open'
            clickType = get(gcf,'SelectionType');
            
            % Check which pool we are looking atr
            if(strcmp(region,'V')),
                responseTrace = activity.V_firing_history(neuron, :, period, epoch);
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

                figure;
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
                        rectangle('Position', [timeToTimeStep(lastStimOnsetTime,dt), 0.001, timeToTimeStep(stimDuration,dt),  1],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);

                        lastStimOnsetTime = [];
                    end
                end

                if(~isempty(lastStimOnsetTime))
                    rectangle('Position', [timeToTimeStep(lastStimOnsetTime, dt), 0.001, timeToTimeStep(numTimeSteps-lastStimOnsetTime, dt), 1],'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]);
                end

                plot(0:(numTimeSteps-1), responseTrace, 'b');
                if ~isempty(s), plot([s s],[0 1],'r'); end % Saccade times
                hXLabel = xlabel('Time (s)');
                hYLabel = ylabel('Firing Rate');
                ylim([0 1.1]);

                xTick = 11:10:numTimeSteps;

                for i=1:length(xTick),
                    xTickLabels{i} = [num2str(stepToTime(xTick(i), dt))];
                end

                set(gca,'XTick', xTick);
                set(gca,'XTickLabel', xTickLabels);
                set(gca,'YTick', [0 1]);

                set([hYLabel hXLabel], 'FontSize', 20);
                set([gca], 'FontSize', 18);
                daspect([40 1 1]);

                box on;
                
            else % weight vectors!
                
                if(strcmp(region,'V')),
                    
                    V_to_C_weights = network.V_to_C_weights;
                    
                    maxSynapseWeight = max(max(V_to_C_weights));
                    
                    S = CLayerProbleFileAnalysis.CLabeProbe_Neurons_S;
                    V = CLayerProbleFileAnalysis.CLabeProbe_Neurons_V;
                    
                    min_S = min(S);
                    max_S = max(S);
                    
                    min_V = min(V);
                    max_V = max(V);
                    
                    min_ = min(min_S, min_V);
                    max_ = max(max_S, max_V);
                    
                    numCSynapses = length(S);

                    figure;
                    hold on;

                    for i=1:numCSynapses,
                        synapticWeights = V_to_C_weights(i,neuron);
                        plot(V(i),S(i),'o','Color',[1, 1 -  synapticWeights/maxSynapseWeight, 1 -  synapticWeights/maxSynapseWeight]);
                    end
                    
                    xlim([min_ max_]);
                    ylim([min_ max_]);
                    
                    hXLabel = xlabel('Retinal Location (deg)');
                    hYLabel = ylabel('Saccade Target (deg)');
                    set([hYLabel hXLabel], 'FontSize', 20);
                    set([gca], 'FontSize', 18);
                   
                %elseif(strcmp(region,'R')),
                %    
                elseif(strcmp(region,'S')),
                    %responseTrace = activity.S_firing_history(neuron, :, period, epoch);
                elseif(strcmp(region,'C')),
                    
                    %% START plotting!
                    
                    C_to_R_weights = network.C_to_R_weights;
                    
                    weightvector = C_to_R_weights(:, neuron);
                    
                    figure;
                    
                    plot(weightvector);
                end

            end
        end

    end


end

