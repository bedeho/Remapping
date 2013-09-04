
dt = 0.005;

Testing_StimuliControl('basic', dt);
Testing_SaccadeControl('basic', dt);
Testing_Kusonoki('basic', dt);
Testing_DuhamelTruncation('basic', dt);
Testing_DuhamelRemappingTrace('basic', dt);
Testing_DuhamelRemapping('basic', dt);
Testing_CLayerProbeTask('basic', dt);

GenerateExperiment(dt);