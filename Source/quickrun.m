    
    %% Import global variables
    declareGlobalVars();
    
    global STIMULI_FOLDER;
    global EXPERIMENTS_FOLDER;
    
    experimentName = 'noCRF-sparse0.2';
    
    experimentFolder = [EXPERIMENTS_FOLDER experimentName filesep];
    
    %% DO NOT CHANGE, basic2 has fixed dt
    %dt = 0.010; 
    dt = 0.005;
    %dt = 0.002;
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
            Testing_StimuliControl('basic2', dt, Training_RF_Locations+Training_Saccades);
            Testing_SaccadeControl('basic2', dt, Training_Saccades);
        %end

        Testing_CLayerProbeTask('basic', dt);

        Testing_DuhamelRemapping('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelTruncation('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_DuhamelRemappingTrace('basic', dt, Training_RF_Locations, Training_Saccades);
        
        Testing_Kusonoki('basic', dt, Training_RF_Locations, Training_Saccades);

    end
    
    %stimulinames = {'basic-StimuliControl','basic-CLayerProbe'};

    %{
    stimulinames = {'basic-StimuliControl', ...
                    'basic-SaccadeControl', ...
                    'basic2-StimuliControl', ...
                    'basic2-SaccadeControl', ...
                    'basic-DuhamelRemapping', ...
                    'basic-DuhamelRemappingTrace', ... 
                    'basic-DuhamelTruncation', ...
                    'basic-CLayerProbe', ...
                    'basic-Kusonoki'};
    %}
    
        stimulinames = {'basic-StimuliControl', ...
                    'basic-SaccadeControl', ...
                    'basic2-StimuliControl', ...
                    'basic2-SaccadeControl', ...
                    'basic-DuhamelRemapping', ...
                    'basic-DuhamelRemappingTrace', ... 
                    'basic-DuhamelTruncation'
                    };
    
    %'basic2-StimuliControl', ... % dt=0.01
    %'basic2-SaccadeControl', ... % dt=0.01      

    %% Run

    GenerateExperiment(experimentName, dt, stimulinames, trainingStimuli);
    
    %% Backup 

    % Compress source code folder
    system(['tar -cjvf ' experimentFolder 'source.tbz .']);
    
    % Iterate stimuli and move to experiment
    for s=1:length(stimulinames),
        
        % Get stimuli name
        stimName = stimulinames{s};
        
        % Get stimuli folder
        stimuliFolder = [STIMULI_FOLDER stimName];
        
        % Compress stimuli folder
        system(['tar -cjvPf ' experimentFolder stimName '.tbz ' stimuliFolder]);
    
    end