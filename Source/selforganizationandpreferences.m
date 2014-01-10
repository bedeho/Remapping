
%
%  selforganizationandpreferences.m
%  Remapping
%
%  Created by Bedeho Mender on 09/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function selforganizationandpreferences()

    % Blank
    CLayerProbleFile_blank = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/baseline-onsettune/S_presaccadic_onset=0.07/BlankNetwork/analysis-basic-CLayerProbe.mat';
    CLayerProbleFileAnalysis_blank = load(CLayerProbleFile_blank);
    S_blank = CLayerProbleFileAnalysis_blank.CLabeProbe_Neurons_S;
    V_blank = CLayerProbleFileAnalysis_blank.CLabeProbe_Neurons_V;
    
    % Trained
    CLayerProbleFile_trained = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/baseline-onsettune/S_presaccadic_onset=0.07/TrainedNetwork/analysis-basic-CLayerProbe.mat';
    CLayerProbleFileAnalysis_trained = load(CLayerProbleFile_trained);
    S_trained = CLayerProbleFileAnalysis_trained.CLabeProbe_Neurons_S;
    V_trained = CLayerProbleFileAnalysis_trained.CLabeProbe_Neurons_V;

    % S
    
    [receptivefieldPlot, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms({S_blank}, {S_trained}, 'XTitle', 'Untrained saccade preference (deg)', 'YTitle', 'Trained saccade preference (deg)');
    
    length(S_blank)
    remove_from_S = isnan(S_blank) | isnan(S_trained);
    S_blank(remove_from_S) = [];
    S_trained(remove_from_S) = [];
    length(S_blank)
    
    corrcoef(S_blank,S_trained)
    
    % R
    
    [receptivefieldPlot, yProjectionAxis, scatterAxis, xProjectionAxis, XLim, YLim] = scatterPlotWithMarginalHistograms({V_blank}, {V_trained}, 'XTitle', 'Untrained visual preference (deg)', 'YTitle', 'Trained visual preference (deg)');
    
    length(V_blank)
    remove_from_V = isnan(V_blank) | isnan(V_trained);
    V_blank(remove_from_V) = [];
    V_trained(remove_from_V) = [];
    length(V_blank)
    
    corrcoef(V_blank,V_trained)
end

