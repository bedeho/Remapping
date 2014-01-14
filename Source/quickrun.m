
function quickrun(DELAY)

    %experimentName = ['baseline-delay' num2str(DELAY) ];
    %for d=[0.050 0.100 0.200 0.250], quickrun(d); end;
    
    experimentName = 'LHeiser';
    
    %experimentName = 'prewired';
    %experimentName = 'baseline';
    %experimentName = 'baseline_denser';
    
    
    %% DO NOT CHANGE, basic2 has fixed dt
    %dt = 0.010; 
    %dt = 0.005;
    %dt = 0.004;
    dt = 0.002;
    %dt = 0.001;
    
    %% Setup stimuli
    
    if true,

        % Training
        %Training_Basic('basic', dt);
        [Training_RF_Locations, Training_Saccades, trainingStimuli] = Training_Coordinated('basic', dt);
        [LHeiser_Training_RF_Locations, LHeiser_trainingStimuli] = Training_LHeiser('basic', dt);

        % Testing
        Testing_StimuliControl('basic', dt);
        Testing_SaccadeControl('basic', dt);
        
        %if length(Training_RF_Locations) == 1,
        %    Testing_StimuliControl('basic2', dt, Training_RF_Locations+Training_Saccades);
        %    Testing_SaccadeControl('basic2', dt, Training_Saccades);
        %end

        Testing_CLayerProbeTask('basic', dt);

        % Cont: Remapping
        Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % LHeiser
        Testing_DuhamelRemapping('basic', dt, LHeiser_Training_RF_Locations, LHeiser_trainingStimuli, 'LHeiser', 0.600, 0.100, 0.100, 0.300);
        
        % Truncation
        Testing_DuhamelTruncation('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % Trace
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
        % Kusonoki
        Testing_Kusonoki('basic', dt, Training_RF_Locations, Training_Saccades);
        
    end
    
    %stimulinames = {'basic-StimuliControl','basic-CLayerProbe'};
    
    %'basic2-StimuliControl', ...
    %'basic2-SaccadeControl', ...

    
    stimulinames = {'basic-StimuliControl', ...
                    'basic-SaccadeControl', ...
                    'basic-CLayerProbe', ...
                    'basic-DuhamelRemapping', ...
                    'basic-DuhamelRemappingTrace', ... 
                    'basic-DuhamelTruncation', ...
                    'basic-Kusonoki', ...
                    'basic-LHeiser'}; % , ... 'basic-Kusonoki'
   
    
    %% Run
    %GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli, DELAY);
    
    %GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli);
    GenerateExperiment(experimentName, dt, stimulinames, LHeiser_trainingStimuli);

end