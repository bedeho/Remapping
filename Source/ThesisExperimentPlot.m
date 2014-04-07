
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
        
    %% Prewired
    %{
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
        %rem = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemapping.mat']);
        %remapping_results{i} = rem.DuhamelRemapping_Result;
        
        % Load Duhamel Trace file
        duhamelRemapping = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemappingTrace.mat']);
        remappingtrace_results{i} = duhamelRemapping.DuhamelRemappingTrace_Result;

    end
    %}
    
    %% Selforganizing
    
    experiment  = 'baseline';

    simulationFolder{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork'];
    Legends{1}          = 'Untrained';
    FaceColors{1}       = [67,82,163]/255;
    
    simulationFolder{2} = [EXPERIMENTS_FOLDER experiment '/baseline/TrainedNetwork'];
    Legends{2}          = 'Trained';
    FaceColors{2}       = [238,48,44]/255;
    
    
    %% baseline-delay
    %{
    %experiment  = 'baseline-delay0.3';
    experiment  = 'prewired-delay0.05';

    simulationFolder{1} = [EXPERIMENTS_FOLDER experiment '/baseline/BlankNetwork'];
    Legends{1}          = 'Random'; %'Untrained';
    FaceColors{1}       = [67,82,163]/255;
    
    simulationFolder{2} = [EXPERIMENTS_FOLDER experiment '/baseline/PrewiredNetwork']; %TrainedNetwork
    Legends{2}          = 'Manual'; %'Trained';
    FaceColors{2}       = [238,48,44]/255;
    %}
        
    % Iterate Simulations
    for i=1:length(simulationFolder),
        
        % Load Remapping file
        %rem = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemapping.mat']);
        %remapping_results{i} = rem.DuhamelRemapping_Result;
        
        % Load Duhamel Trace file
        duhamelRemapping = load([simulationFolder{i} filesep 'analysis-basic-DuhamelRemappingTrace.mat']);
        remappingtrace_results{i} = duhamelRemapping.DuhamelRemappingTrace_Result;

    end
    
    
    % Remapping
    %disp('Continous =======');
    %[remLatFig, remScatFig, indexFig] = dump(remapping_results);
    
    % Trace
    disp('Flashed ===========');
    [remLatFig, remScatFig, indexFig] = dump(remappingtrace_results);

    function [f1, f2, f3] = dump(res)
        
        [f1, f2, f3] = remappingPlots(res, FaceColors, Legends);

        for j=1:length(simulationFolder),
            
            latencies = [res{j}.remappingLatency];
            NonNaNLatencies = latencies;
            NonNaNLatencies(isnan(NonNaNLatencies)) = [];
            
            Num_NonNaNLatencies = length(NonNaNLatencies);
            TOTAL = length(latencies);
            
            disp(Legends{j});
            disp(['Valid latencies: ' num2str(Num_NonNaNLatencies) '/' num2str(TOTAL) ]);
            disp(['Average Remapping Latency: ' num2str(mean(1000*NonNaNLatencies))]);
            disp(['Average Remapping Index: ' num2str(mean([res{j}.remapping_index]))]);
            
            x=nnz([res{j}.stimLatency] > [res{j}.remappingLatency]);
            disp(['Predictive: ' num2str(x) '=' num2str(100*x/TOTAL) '%']);
            
            x = nnz([res{j}.remappingLatency] < 0);
            disp(['Presaccadic: ' num2str(x) '=' num2str(100*x/TOTAL) '%']);
            
            disp(['MAX Index: ' num2str(max([res{j}.remapping_index]))]);
            disp(['Stim Latency: ' num2str(mean([res{j}.stimLatency]))]);
            disp(' ');
        end
        
    end
    
end