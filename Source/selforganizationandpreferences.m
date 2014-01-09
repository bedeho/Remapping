
%
%  selforganizationandpreferences.m
%  Remapping
%
%  Created by Bedeho Mender on 09/01/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

function selforganizationandpreferences()

    % Blank
    CLayerProbleFile_blank = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/baseline-onsettune/S_presaccadic_onset=0.07/BlankNetwork/analysis-basic-CLayerProbe.mat'
    CLayerProbleFileAnalysis_blank = load(CLayerProbleFile_blank);
    S_blank = CLayerProbleFileAnalysis_blank.CLabeProbe_Neurons_S;
    V_blank = CLayerProbleFileAnalysis_blank.CLabeProbe_Neurons_V;

    % Trained
    CLayerProbleFile_trained = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/Experiments/baseline-onsettune/S_presaccadic_onset=0.07/TrainedNetwork/analysis-basic-CLayerProbe.mat'
    CLayerProbleFileAnalysis_trained = load(CLayerProbleFile_trained);
    S_trained = CLayerProbleFileAnalysis_trained.CLabeProbe_Neurons_S;
    V_trained = CLayerProbleFileAnalysis_trained.CLabeProbe_Neurons_V;
    
    R_eccentricity = 45;
    S_eccentricity = 30;
    
    % Plot
    plot(S_blank,S_trained,'o');
    
    


end

