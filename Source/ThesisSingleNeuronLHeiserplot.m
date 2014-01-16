
%
%  ThesisSingleNeuronLHeiserplot.m
%  Remapping
%
%  Created by Bedeho Mender on 16/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function ThesisSingleNeuronLHeiserplot()

    % Import global variables
    declareGlobalVars();
    global EXPERIMENTS_FOLDER;
        
    simulation = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.1/TrainedNetwork/analysis-basic-LHeiser.mat'];
    stimuli = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/STIM-basic-LHeiser'];
    
    %% Load stimuli
    sim = load(simulation);
    stim = load([stimuli '/stim.mat']);
    
    
    
    
    %{
    numberOfDirections = 4;
    numberOfUnique = 17;
    xTitle = '\phi^\text{C}';
    
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
        figure('Units','pixels','position', [1000 1000 400 200]); % [1000 1000 400 400]

        h = hist(uniqueResponseCount,0:numberOfDirections)/numberOfUnique;
        bar(0:numberOfDirections, h, 0.7);
        ylim([0 1.01]);
        
        hXLabel = xlabel('Number of saccade directions with remapping');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'YTick', [0 1]);
        pbaspect([0.6 0.3 1])
        %axis square;
        
        % Save plot for bar plot
        results(i,:) = h;
        
    end

    %}
end