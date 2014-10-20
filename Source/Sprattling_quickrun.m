
function Sprattling_quickrun()

    %for d=[0.050 0.100 0.150 0.200 0.250 0.300],
    %for d=0.05:0.1:1,
    
    % DILUTION
    % d=0.05
        
        %experimentName = ['baseline-delay' num2str(d) ];
        %experimentName = ['prewired-delay' num2str(d) ];
        
        %experimentName = 'V_dilution';
        %experimentName = 'S_dilution';
        %experimentName = ['dilution-' num2str(d)];
        
        %experimentName = 'S_onsetvariability';
        
        %experimentName = 'prewired';
        
        % sprattlig
        %experimentName = 'sprattling_nonplastic';
        experimentName = 'sprattling_visual_learning_bigepoch20-connectivitfix-tune53'; %tune40
        
        %% DO NOT CHANGE, basic2 has fixed dt
        %dt = 0.010;
        %dt = 0.005;
        %dt = 0.004;
        dt = 0.002;
        %dt = 0.001;
        
        %% Setup stimuli
        
        % Training
        %Training_Basic('basic', dt);
        [Training_RF_Locations, Training_Saccades, trainingStimuli] = Training_Coordinated('basic', dt);
        
        % Testing
        Testing_StimuliControl('basic', dt);
        Sprattling_Testing_StimuliControl('basic', dt);
        Testing_SaccadeControl('basic', dt);
        
        %if length(Training_RF_Locations) == 1,
        %    Testing_StimuliControl('basic2', dt, Training_RF_Locations+Training_Saccades);
        %    Testing_SaccadeControl('basic2', dt, Training_Saccades);
        %end
        
        Testing_CLayerProbeTask('basic', dt);
        
        % Cont: Remapping
        %Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
        
        stimulinames = {'basic-StimuliControl', ...
            'basic-SaccadeControl', ...
            'basic-DuhamelRemappingTrace', ...
            'basic-Sprattling_StimuliControl'
            };
        
        %'basic-CLayerProbe', ...
        
        %% Run
        Sprattling_GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli);
    %end
   
end