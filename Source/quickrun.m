
    %% DO NOT CHANGE, basic2 has fixed dt
    %dt = 0.010; 
    dt = 0.005;

    if true,

        %% Training
        %Training_Basic('basic', dt);

        [Training_RF_Locations, Training_Saccades, trainingStimuli] = Training_Coordinated('basic', dt)

        %% Testing

        %Stim/Sacc controls
        Testing_StimuliControl('basic', dt);
        Testing_SaccadeControl('basic', dt);

        %CLayerProbe task
        Testing_CLayerProbeTask('basic', dt);

        % Perisaccadic testing tasks

        Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelTruncation('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_Kusonoki('basic', dt, Training_RF_Locations, Training_Saccades);

    end
    
    stimulinames = {'basic-StimuliControl', ...
                    'basic2-StimuliControl', ... % dt=0.01
                    'basic-SaccadeControl', ...
                    'basic2-SaccadeControl', ... % dt=0.01
                    'basic-DuhamelRemapping', ...
                    'basic-DuhamelRemappingTrace', ... 
                    'basic-DuhamelTruncation', ...
                    'basic-CLayerProbe', ...
                    'basic-Kusonoki'};

    GenerateExperiment('learning16',dt, stimulinames, trainingStimuli);