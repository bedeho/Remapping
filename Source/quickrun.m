
function quickrun()

    experimentName = 'prewired'; % baseline-multiRF3-tuneC-to-Rpsi
    
    
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

        % Testing
        Testing_StimuliControl('basic', dt);
        Testing_SaccadeControl('basic', dt);
        
        %if length(Training_RF_Locations) == 1,
        %    Testing_StimuliControl('basic2', dt, Training_RF_Locations+Training_Saccades);
        %    Testing_SaccadeControl('basic2', dt, Training_Saccades);
        %end

        Testing_CLayerProbeTask('basic', dt);

        Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelTruncation('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
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
                    'basic-Kusonoki'}; % , ... 'basic-Kusonoki'
    
    %% Run
    GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli);

end