
function quickrun_LHeiser()
    
    experimentName = 'LHeiser';

    %% DO NOT CHANGE, basic2 has fixed dt
    %dt = 0.010; 
    %dt = 0.005;
    %dt = 0.004;
    dt = 0.002;
    %dt = 0.001;
    
    %% Setup stimuli

    % Training
    [LHeiser_Training_RF_Locations, LHeiser_Saccades, LHeiser_trainingStimuli] = Training_LHeiser('basic', dt);

    % Testing
    Testing_StimuliControl('basic', dt);
    Testing_SaccadeControl('basic', dt);
    Testing_CLayerProbeTask('basic', dt);

    % LHeiser
    Testing_DuhamelRemapping('basic', dt, LHeiser_Training_RF_Locations, LHeiser_Saccades, 'LHeiser', 0.600, 0.100, 0.100, 0.300);

    %% Run
    GenerateExperiment(experimentName, dt, {'basic-StimuliControl', 'basic-SaccadeControl', 'basic-CLayerProbe', 'basic-LHeiser'}, LHeiser_trainingStimuli);

end