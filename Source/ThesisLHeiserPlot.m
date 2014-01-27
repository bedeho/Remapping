
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
    %{
    simulationFolder{1} = [EXPERIMENTS_FOLDER 'LHeiser/baseline'];
    Legends{1} = 'value 1';
    
    %simulationFolder{2} = [EXPERIMENTS_FOLDER 'LHeiser/baseline'];
    %Legends{2} = 'value 2';
    
    numberOfDirections = 3;
    numberOfUnique = 2;
    xTitle = 'The x tite';
    %}
    
    %% LHeiser_C_to_R_connectivity
    
    simulationFolder{1} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.1'];
    Legends{1}          = '0.1';
    simulationFolder{2} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.2'];
    Legends{2}          = '0.2';
    simulationFolder{3} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.3'];
    Legends{3}          = '0.3';
    simulationFolder{4} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.4'];
    Legends{4}          = '0.4';
    simulationFolder{5} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.5'];
    Legends{5}          = '0.5';
    simulationFolder{6} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.6'];
    Legends{6}          = '0.6';
    simulationFolder{7} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.7'];
    Legends{7}          = '0.7';
    simulationFolder{8} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.8'];
    Legends{8}          = '0.8';
    simulationFolder{9} = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=0.9'];
    Legends{9}          = '0.9';
    simulationFolder{10}= [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=1'];
    Legends{10}         = '1.0';
    
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

        %% Number of locations responsiveness distribution
        
        %{
        
        % Compute        
        diffUniqueIndexes = trained_uniqueIndexes - blank_uniqueIndexes;
        uniqueResponseCount = sum(diffUniqueIndexes > 0);
        
        % Show figure
        figure('Units','pixels','position', [1000 1000 400 200]); % [1000 1000 400 400]

        h = hist(uniqueResponseCount,0:numberOfDirections)/numberOfUnique;
        bar(0:numberOfDirections, h, 0.7);
        ylim([0 1.01]);
        
        hXLabel = xlabel('Number of locations with remapping');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'YTick', [0 1]);
        pbaspect([0.6 0.3 1])
        %axis square;
        
        % Save plot for bar plot
        results(i,:) = h;
        
        
        %}
        
        %% selectivity distribution
        %{
        % Make histograms
        SI_untrained = [blank_network.LHeiserAnalysis.SI];
        SI_trained = [trained_network.LHeiserAnalysis.SI];

        X = 0:0.05:1;

        %h_untrained = histc(SI_untrained, X);
        h_trained = histc(SI_trained, X);

        H = [ h_trained(1:end)']; % h_untrained(1:end)'

        figure('Units','pixels','position', [1000 1000 400 200]);

        hBar = bar(X,H,1.0);  %bar(X,H,1.0,'stacked','LineStyle','none'); 
        set(hBar(1),'FaceColor', 'r');
        %set(hBar(2),'FaceColor', 'r');

        xlim([-0.1 1.1])

        hXLabel = xlabel('Selectivity Index');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
            %}
        
        %% preference distribution
        %{
        % Make histograms
        pref_trained = [trained_network.LHeiserAnalysis.retinal_preference];

        X = -45:5:45;

        h_trained = histc(pref_trained, X);

        figure('Units','pixels','position', [1000 1000 400 200]);

        hBar = bar(X,h_trained,1.0); 

        xlim([-50 50])

        hXLabel = xlabel('Preferred Location (deg)');
        hYLabel = ylabel('Frequency');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'XTick',-45:10:45);
        %}
        
        %% Number of locations vs. max remapping index
        %{
        % Compute        
        diffUniqueIndexes = trained_uniqueIndexes - blank_uniqueIndexes;
        uniqueResponseCount = sum(diffUniqueIndexes > 0);
        maxIndexPerNeuron = max(trained_uniqueIndexes);
        
        means = zeros(1,numberOfDirections);
        err = zeros(1,numberOfDirections);
        
        for j=1:numberOfDirections,
            
            v = maxIndexPerNeuron(uniqueResponseCount == j);
            means(j) = mean(v);
            err(j) = std(v);
        end
        
        % Show figure
        figure('Units','pixels','position', [1000 1000 400 200]); % [1000 1000 400 400]
        
        errorbar(1:numberOfDirections, means, err);
        
        ylim([0 sqrt(2)]);
        xlim([0 (numberOfDirections+1)]);
        
        hXLabel = xlabel('Number of locations with remapping');
        hYLabel = ylabel('Maximum Remapping Index');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'YTick', [0 1], 'XTick', 1:numberOfDirections);
        pbaspect([0.6 0.3 1])
        %}
        
        %% RI,SI as function of RF
        
        meanRI = mean(trained_uniqueIndexes);
        errorRI = std(trained_uniqueIndexes);
        RF = trained_network.LHeiserAnalysis.Unique_Training_RF_Locations;
        SI = trained_network.LHeiserAnalysis.SI;

        % Show figure
        figure('Units','pixels','position', [1000 1000 400 200]); % [1000 1000 400 400]
        
        errorbar(RF, meanRI, errorRI, 'o');
        hold on;
        plot(RF,SI,'*g');
        
        YLim = [0 sqrt(2)];
        XLim = [-45 45];
        
        ylim(YLim);
        xlim(XLim);
        
        hXLabel = xlabel('Retinal Location (deg)');
        hYLabel = ylabel('');
        
        if(i==1),
            legend('Mean Remapping Index','Selectivity Index');
        end
        
        %legend('boxoff');
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        set(gca, 'YTick', [0 1], 'XTick', -45:10:45);
        pbaspect([0.6 0.3 1])
                
    end
    
    %{
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
    %}
    
end