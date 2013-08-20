
%
%  viewNeuronDynamics.m
%  Remapping
%
%  Created by Bedeho Mender on 26/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function viewNeuronDynamics(activityFile, stimuliName)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;
    
    if nargin == 0,
        activityFile    = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/prewired/baseline/PrewiredNetwork/activity-basic-KusonokiTesting.mat';
        stimuliName     = 'basic-KusonokiTesting';
    end
    
    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load([STIMULI_FOLDER stimuliName filesep 'stim.mat']);
    
    % Set parameters
    R_N = activity.R_N;
    C_N = activity.C_N;
    dt = stimuli.dt;
    
    % Make figure
    figure('name',activityFile,'Position', [100, 100, 1049, 895]);
    
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
        
        % If C is empty, just fill with blank
        if ~isempty(activity.C_firing_history),
            C_firingrate = activity.C_firing_history(:, :, period, epoch);
            C_activation = activity.C_activation_history(:, :, period, epoch);
            
            includeC_layer = 1;
            numRows = 2+5;
            
        else
            C_firingrate = zeros(C_N, activity.numPeriods);
            C_activation = zeros(C_N, activity.numPeriods);
            
            includeC_layer = 0;
            numRows = 5;
        end
        
        % Plot
        s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);
        numTimeSteps = length(stimuli.stimuli{period}.eyePositionTrace);

        subplot(6,2,1);
        imagesc(E_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('E Firing');

        %{
        subplot(numRows,2,2);
        imagesc(E_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('E Firing');
        %}
        
        subplot(numRows,2,3);
        imagesc(V_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('V Firing');

        %{
        subplot(numRows,2,4);
        imagesc(V_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('V Firing');
        %}
        
        subplot(numRows,2,5);
        imagesc(R_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Firing');

        subplot(numRows,2,6);
        imagesc(R_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Activation');

        subplot(numRows,2,7);
        imagesc(S_firingrate);
        colorbar
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('S Firing');

        subplot(numRows,2,8);
        imagesc(S_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('S Activation');

        if(includeC_layer),
            subplot(numRows,2,[9 11]);
            imagesc(C_firingrate);
            hold on;colorbar
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Firing');

            subplot(numRows,2,[10 12]);
            imagesc(C_activation);
            hold on;colorbar
            if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
            colorbar
            title('C Actiation');
            
            nextplot = 13;
        else
            nextplot = 9;
        end
        
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
        
    end
end

