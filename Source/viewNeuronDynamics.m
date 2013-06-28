
%
%  viewNeuronDynamics.m
%  Remapping
%
%  Created by Bedeho Mender on 26/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function viewNeuronDynamics(activityFile, stimuliFile)

    % Load input files
    disp('Loading input files...');
    activity = load(activityFile);
    stimuli  = load(stimuliFile);
    
    % Set parameters
    R_N = activity.R_N;
    S_N = activity.S_N;
    C_N = activity.C_N;
    eyePositionTrace = stimuli.eyePositionTrace;
    retinalTargetTraces = stimuli.retinalTargetTraces;
    
    % Visualization loop
    figure('Position', [100, 100, 1049, 895]);
    
    while(true),
        
        % Get desired epoch or quit
        if activity.numEpochs > 1,
            str = input(['Epoch # (1, ' num2str(activity.numEpochs) '), or 0 to quit'],'s');
            epoch = num2str(str);
            
            % Quit if we are done
            if epoch == 0,
                return;
            end
        end
        
        % Get desired period
        str = input(['Period # (1, ' num2str(activity.numPeriods) ')'],'s');
        period = num2str(str);
        
        % Load
        V_firingrate = activity(epoch, period, 1).V_firingrate;
        R_firingrate = activity(epoch, period, 1).R_firingrate;
        S_firingrate = activity(epoch, period, 1).S_firingrate;
        C_firingrate = activity(epoch, period, 1).C_firingrate;

        V_activation = activity(epoch, period, 2).V_activation;
        R_activation = activity(epoch, period, 2).R_activation;
        S_activation = activity(epoch, period, 2).S_activation;
        C_activation = activity(epoch, period, 2).C_activation;
        
        % Plot
        s = timeToTimeStep(stimuli{period}.saccadeTimes);

        subplot(5,2,1);
        imagesc(flipud(V_firingrate));
        hold on;plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r');
        colorbar
        title('V Firing');

        subplot(5,2,2);
        imagesc(flipud(V_activation));
        hold on;plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r');
        colorbar
        title('V Firing');

        subplot(5,2,3);
        imagesc(flipud(R_firingrate));
        hold on;plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r');
        colorbar
        title('R Firing');

        subplot(5,2,4);
        imagesc(flipud(R_activation));
        hold on;plot([s s],[ones(R_N,1) R_N*ones(R_N,1)],'r');
        colorbar
        title('R Activation');

        subplot(5,2,5);
        imagesc(flipud(S_firingrate));
        colorbar
        hold on;plot([s s],[ones(S_N,1) S_N*ones(S_N,1)],'r');
        title('S Firing');

        subplot(5,2,6);
        imagesc(flipud(S_activation));
        hold on;plot([s s],[ones(S_N,1) S_N*ones(S_N,1)],'r');
        colorbar
        title('S Activation');

        subplot(5,2,7);
        imagesc(flipud(C_firingrate));
        hold on;plot([s s],[ones(C_N,1) C_N*ones(C_N,1)],'r');
        colorbar
        title('C Firing');

        subplot(5,2,8);
        imagesc(flipud(C_activation));
        hold on;plot([s s],[ones(C_N,1) C_N*ones(C_N,1)],'r');
        colorbar
        title('C Actiation');

        subplot(5,2,9);
        plot(eyePositionTrace, 'r');
        hold on;
        plot(retinalTargetTraces' , 'b');
        xlabel('Time step');
        legend({'Eye Position','Stimuli Retinal Locations'})

        subplot(5,2,10);
        plot(eyePositionTrace, 'r');
        hold on;
        plot(retinalTargetTraces' , 'b');
        xlabel('Time step');
        legend({'Eye Position','Stimuli Retinal Locations'})    

        % ylim(min(min(eyePositionTrace),min(min(retinalTargetTraces)) max(max(eyePositionTrace),max(max(retinalTargetTraces)))]);
    end
end