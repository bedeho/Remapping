
%
%  ThesisLHeiserPlot.m
%  Remapping
%
%  Created by Bedeho Mender on 01/15/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function ThesisLHeiserPlot()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
        
    % Simulations 
    simulationFolder{1} = [EXPERIMENTS_FOLDER 'LHeiser/baseline'];
    Legends{1} = 'value 1';
    
    %simulationFolder{2} = [EXPERIMENTS_FOLDER 'LHeiser/baseline'];
    %Legends{2} = 'value 2';
    
    numberOfDirections = 3;
    numberOfUnique = 2;
    xTitle = 'The x tite';
    
    % Iterate simulations
    numSimulations = length(simulationFolder);
    results = zeros(numSimulations, numberOfDirections + 1);
    
    % Iterate Simulations
    for i=1:numSimulations,
        
        % Load untrained
        blank_network = load([simulationFolder{i} filesep 'BlankNetwork/analysis-basic-LHeiser.mat']);
        blank_uniqueIndexes = blank_network.LHeiserAnalysis.uniqueIndexes
        
        % Load trained
        trained_network = load([simulationFolder{i} filesep 'TrainedNetwork/analysis-basic-LHeiser.mat']);
        trained_uniqueIndexes = trained_network.LHeiserAnalysis.uniqueIndexes

        % Compute        
        diffUniqueIndexes = trained_uniqueIndexes - blank_uniqueIndexes;
        uniqueResponseCount = sum(diffUniqueIndexes > 0);
        
        % Show figure
        figure;

        h = hist(uniqueResponseCount,0:numberOfDirections)/numberOfUnique;
        bar(0:numberOfDirections, h, 0.7);
        ylim([0 1.01]);
        
        hXLabel = xlabel('Number of saccade directions with remapping');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'YTick', [0 1]);
        
        axis square;
        
        % Save plot for bar plot
        results(i,:) = h;
        
    end
    
    % Plot population plot
    figure;
    plot(results');
    ylim([0 1.01]);
    hXLabel = xlabel(xTitle);
    hYLabel = ylabel('Frequency');
    set([hYLabel hXLabel], 'FontSize', 14);
    set(gca, 'FontSize', 12);
    set(gca, 'YTick', [0 1]);
    
    legend(Legends);
    
end