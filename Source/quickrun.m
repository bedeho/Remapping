
%% DO NOT CHANGE, basic2 has fixed dt
%dt = 0.010; 
dt = 0.005;

if true,
    Training_Basic('basic', dt);

    Testing_StimuliControl('basic', dt);
    Testing_SaccadeControl('basic', dt);
    Testing_Kusonoki('basic', dt);
    Testing_DuhamelTruncation('basic', dt);
    Testing_DuhamelRemappingTrace('basic', dt);
    Testing_DuhamelRemapping('basic', dt);
    Testing_CLayerProbeTask('basic', dt);
end

GenerateExperiment('learning13',dt);