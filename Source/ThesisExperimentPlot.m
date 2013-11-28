
%
%  ThesisExperimentPlot.m
%  Remapping
%
%  Created by Bedeho Mender on 25/11/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function ThesisExperimentPlot()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
    
    % Basic
    experiment  = 'test_plotting';
    
    simulationFolder{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork/'];
    Legends{1}          = 'Untrained';
    FaceColors{1}       = [67,82,163]/255;
    
    simulationFolder{2} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/'];
    Legends{2}          = 'Trained';
    FaceColors{2}       = [238,48,44]/255;
    
    % Iterate Simulations
    for i=1:length(simulationFolder),
        
        % Load Remapping file
        rem = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemapping.mat']);
        remapping_results{i} = rem.DuhamelRemapping_Result;
        
        % Load Duhamel Trace file
        duhamelRemapping = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemappingTrace.mat']);
        remappingtrace_results{i} = duhamelRemapping.DuhamelRemappingTrace_Result;

    end
    
    % Perform plots
    [remLatFig, remScatFig, indexFig] = remappingPlots(remapping_results, FaceColors, Legends);
    [remLatFig, remScatFig, indexFig] = remappingPlots(remappingtrace_results, FaceColors, Legends);
    
end

%{ 
old check existance BS
if(exist(remFile, 'file')),

    rem = load(remFile);
    remapping_results{i} = rem.DuhamelRemapping_Result;

end
%}

