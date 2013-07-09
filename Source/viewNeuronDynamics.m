
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
    
    V_firing_history = activity.V_firing_history;
    R_firing_history = activity.R_firing_history;
    S_firing_history = activity.S_firing_history;
    C_firing_history = activity.C_firing_history;

    V_activation_history = activity.V_activation_history;
    R_activation_history = activity.R_activation_history;
    S_activation_history = activity.S_activation_history;
    C_activation_history = activity.C_activation_history;
    
    % Set parameters
    R_N = activity.R_N;
    C_N = activity.C_N;
    dt = stimuli.dt;
    
    % Make figure
    figure('Position', [100, 100, 1049, 895]);
    
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
        V_firingrate = V_firing_history(:, :, period, epoch);
        R_firingrate = R_firing_history(:, :, period, epoch);
        S_firingrate = S_firing_history(:, :, period, epoch);
        

        V_activation = V_activation_history(:, :, period, epoch);
        R_activation = R_activation_history(:, :, period, epoch);
        S_activation = S_activation_history(:, :, period, epoch);
        
        if ~isempty(C_firing_history),
            C_firingrate = C_firing_history(:, :, period, epoch);
            C_activation = C_activation_history(:, :, period, epoch);
        else
            C_firingrate = zeros(C_N, activity.numPeriods);
            C_activation = zeros(C_N, activity.numPeriods);
        end
        
        % Plot
        s = timeToTimeStep(stimuli.stimuli{period}.saccadeTimes, dt);

        subplot(5,2,1);
        imagesc(V_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('V Firing');

        subplot(5,2,2);
        imagesc(V_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('V Firing');

        subplot(5,2,3);
        imagesc(R_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Firing');

        subplot(5,2,4);
        imagesc(R_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('R Activation');

        subplot(5,2,5);
        imagesc(S_firingrate);
        colorbar
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        title('S Firing');

        subplot(5,2,6);
        imagesc(S_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('S Activation');

        subplot(5,2,7);
        imagesc(C_firingrate);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('C Firing');

        subplot(5,2,8);
        imagesc(C_activation);
        hold on;colorbar
        if ~isempty(s), plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r'); end
        colorbar
        title('C Actiation');

        subplot(5,2,9);
        cla
        plot(stimuli.stimuli{period}.eyePositionTrace, 'r');
        hold on;
        plot(stimuli.stimuli{period}.retinalTargetTraces', 'b');
        xlabel('Time step');
        legend({'Eye Position','Stimuli Retinal Locations'});
        ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
        set(gca,'YDir','reverse');
        
        subplot(5,2,10);
        cla
        plot(stimuli.stimuli{period}.eyePositionTrace, 'r');
        hold on;
        plot(stimuli.stimuli{period}.retinalTargetTraces', 'b');
        xlabel('Time step');
        legend({'Eye Position','Stimuli Retinal Locations'});
        ylim([-45 45]); % we hard code limit since not all stimuli has stimuli.R_eccentricity
        set(gca,'YDir','reverse');
        
    end
end

