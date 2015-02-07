
function quickrun()

    %for d=[0.050 0.100 0.150 0.200 0.250 0.300],
    %for d=0.05:0.1:1,
    
    % DILUTION
    d=0.05
        
        %experimentName = ['baseline-delay' num2str(d) ];
        %experimentName = ['prewired-delay' num2str(d) ];
        
        %experimentName = 'V_dilution';
        %experimentName = 'S_dilution';
        %experimentName = ['dilution-' num2str(d)];
        
        %experimentName = 'S_onsetvariability';
        
        %experimentName = 'prewired';
        
        experimentName = 'classic-lheiser';        
        %experimentName = 'classic';
        
        %% DO NOT CHANGE, basic2 has fixed dt
        %dt = 0.010;
        %dt = 0.005;
        %dt = 0.004;
        dt = 0.002;
        %dt = 0.001;
        
        %% Setup stimuli
        
        % Training
        %[Training_RF_Locations, Training_Saccades, trainingStimuli] = Training_Coordinated('basic', dt);
        
        [Training_RF_Locations, Training_Saccades, trainingStimuli] = Training_LHeiser('basic', dt);
        
        % Testing
        Testing_StimuliControl('basic', dt);
        Testing_SaccadeControl('basic', dt);
        
        
        %%%
        %Sprattling_Testing_StimuliControl('basic', dt);
        %%%
        
        %if length(Training_RF_Locations) == 1,
        %    Testing_StimuliControl('basic2', dt, Training_RF_Locations+Training_Saccades);
        %    Testing_SaccadeControl('basic2', dt, Training_Saccades);
        %end
        
        Testing_CLayerProbeTask('basic', dt);
        
        % Cont: Remapping
        Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % Truncation
        Testing_DuhamelTruncation('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % Trace
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % Kusonoki
        Testing_Kusonoki('basic', dt, Training_RF_Locations, Training_Saccades);

        %% Debug C
        %{
        stimulinames = {'basic-StimuliControl', ...
                        'basic-SaccadeControl'
                        };
        %}
        
        %% classic
        
        stimulinames = {'basic-StimuliControl', ...
                        'basic-SaccadeControl', ...
                        'basic-CLayerProbe', ...
                        'basic-DuhamelRemappingTrace'
                        };
       
                        %, ...
                        %'basic-Kusonoki'

        
        %% Run
        GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli, d);
    %end
   
end