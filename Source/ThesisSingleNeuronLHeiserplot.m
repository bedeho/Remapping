
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
        
    simulation_untrained = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=1/BlankNetwork/analysis-basic-LHeiser.mat'];
    simulation_trained = [EXPERIMENTS_FOLDER 'LHeiser_C_to_R_connectivity/C_to_R_connectivity=1/TrainedNetwork/analysis-basic-LHeiser.mat'];
    
    % Load data
    sim_untrained = load(simulation_untrained);
    sim_trained = load(simulation_trained);
    
    sim_trained.LHeiserAnalysis
    
    % Get variables
    numUnique = sim_trained.LHeiserAnalysis.numUnique;
    retinalSources = sim_trained.LHeiserAnalysis.retinalSources;
    uniqueIndexes_trained = sim_trained.LHeiserAnalysis.uniqueIndexes;
    uniqueIndexes_untrained = sim_untrained.LHeiserAnalysis.uniqueIndexes;
    
    % Plotting variable
    upperYlimit = sqrt(2);
    
    %% Single neuron plots

    % Iterate each neuron
    for i=1:numUnique,
        
        figureHandle = figure('Units','pixels','position', [1000 1000 400 200]); % [1000 1000 400 400]
        
        pref = sim_trained.LHeiserAnalysis.retinal_preference(i);
        
        ret_src = retinalSources(:,i);
        RI_untrained = uniqueIndexes_untrained(:,i);
        RI_trained = uniqueIndexes_trained(:,i);
        
        % sort retinally
        [X,IX] = sort(ret_src);
        RI_untrained = RI_untrained(IX);
        RI_trained = RI_trained(IX);
        
        hold on;
        plot(X, RI_untrained, '-*b'); % untrained
        plot(X, RI_trained, '-or'); % trained

        if(~isnan(pref)),
            plot([pref pref], [-0.01 upperYlimit], '--k'); % preference
        end
        
        text(-20, upperYlimit,['SI = ' num2str(sim_trained.LHeiserAnalysis.SI(i),'%.4f')], 'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14)
        
        xlim([-45 45]);
        ylim([-0.01 upperYlimit]);
        
        hXLabel = xlabel('Retinal location (deg)');
        hYLabel = ylabel('Remapping Index');
        legend('Untrained','Trained');
        legend('boxoff')
        box('on')
        
        set([hYLabel hXLabel], 'FontSize', 14);
        set(gca, 'FontSize', 12);
        
        %% Save
        name = ['C:\Users\Sindre\Desktop\out\LHeiser_RI_' num2str(i)];
        saveas(figureHandle,[name '.eps'], 'epsc');
        saveas(figureHandle,[name '.png']);
        close(figureHandle);
        
    end

end