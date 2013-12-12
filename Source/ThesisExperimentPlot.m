
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
    
    %% Basic
    %{
    experiment  = 'test_plotting';
    
    simulationFolder{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork/'];
    Legends{1}          = 'Untrained';
    FaceColors{1}       = [67,82,163]/255;
    
    simulationFolder{2} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork/'];
    Legends{2}          = 'Trained';
    FaceColors{2}       = [238,48,44]/255;
    %}
    
    %% Prewired
    experiment  = 'prewired';

    simulationFolder{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork/'];
    Legends{1}          = 'Random';
    FaceColors{1}       = [67,82,163]/255;
    
    simulationFolder{2} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork/'];
    Legends{2}          = 'Manual';
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
    
    % Remapping
    disp('Continous =======');
    [remLatFig, remScatFig, indexFig] = dump(remapping_results);
    
    % Trace
    disp('Flashed ===========');
    [remLatFig, remScatFig, indexFig] = dump(remappingtrace_results);

    function [f1, f2, f3] = dump(res)
        
        [f1, f2, f3] = remappingPlots(res, FaceColors, Legends);

        for j=1:length(simulationFolder),
            
            M = length([res{j}.remappingLatency]);
            
            disp(Legends{j});
            disp(['Average Remapping Latency: ' num2str(mean(1000*[res{j}.remappingLatency]))]);
            disp(['Average Remapping Index: ' num2str(mean([res{j}.remapping_index]))]);
            
            x=nnz([res{j}.stimLatency] > [res{j}.remappingLatency]);
            disp(['Predictive: ' num2str(x) '=' num2str(100*x/M) '%']);
            
            x = nnz([res{j}.remappingLatency] < 0);
            disp(['Presaccadic: ' num2str(x) '=' num2str(100*x/M) '%']);
            
            disp(['MAX Index: ' num2str(max([res{j}.remapping_index]))]);
            disp(['Stim Latency: ' num2str(mean([res{j}.stimLatency]))]);
            disp(' ');
        end
        
    end
    
end