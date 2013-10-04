
%% DO NOT CHANGE, basic2 has fixed dt
%dt = 0.010; 
dt = 0.005;

if true,
    
    %% Training
    %Training_Basic('basic', dt);
    
    [Training_RF_Locations, Training_Saccades] = Training_Coordinated('basic', dt)

    %% Testing
    
    %Stim/Sacc controls
    Testing_StimuliControl('basic', dt);
    Testing_SaccadeControl('basic', dt);
    
    %CLayerProbe task
    Testing_CLayerProbeTask('basic', dt);
    
    %
    %Testing_Kusonoki('basic', dt);
    %Testing_DuhamelTruncation('basic', dt);
    %Testing_DuhamelRemappingTrace('basic', dt);
    Testing_DuhamelRemapping('basic', dt);
    
end

GenerateExperiment('learning13',dt);