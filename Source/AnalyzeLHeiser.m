
%
%  AnalyzeLHeiser.m
%  Remapping
%
%  Created by Bedeho Mender on 09/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function [LHeiserAnalysis] = AnalyzeLHeiser(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli)

    % Do underlying analysis
    LHeiser_DuhamelRemappingTrace = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity, stim_stimuli, sacc_control_activity, sacc_stimuli);
    
    % Save indexes for all neurons and directions in 2d-array
    Unique_Training_RF_Locations = unique(stimuli.currentRF);
    numUnique = length(Unique_Training_RF_Locations);
    numberOfDirections = activity.numPeriods/numUnique;
    uniqueIndexes = reshape([LHeiser_DuhamelRemappingTrace.remapping_index], numberOfDirections, numUnique);
    
    % Find retinal source for each trial
    numPeriods = length(stimuli.stimuli);
    
    for i=1:numPeriods,
        r=stimuli.stimuli{i}.headCenteredTargetLocations;
        retinalSources(i) = r;
    end
    
    retinalSources = reshape(retinalSources, numberOfDirections, numUnique);
    
    % normalization term
    sums = sum(uniqueIndexes);
    
    % Location preerence
    retinal_preference = sum(uniqueIndexes.*retinalSources)./sums;
    
    % Selectivity index
    normalized_indexes = uniqueIndexes ./ repmat(sums, numberOfDirections, 1);
    shannon_information_measure = normalized_indexes.*(log(normalized_indexes)/log(2));
    shannon_information_measure(isnan(shannon_information_measure)) = 0;
    entropies = -sum(shannon_information_measure);
    maxEntropy = log(numberOfDirections)/log(2);
    SI = 1-entropies/maxEntropy;

    % Save
    LHeiserAnalysis.numUnique                    = numUnique;
    LHeiserAnalysis.numberOfDirections           = numberOfDirections;
    LHeiserAnalysis.uniqueIndexes                = uniqueIndexes;
    LHeiserAnalysis.Unique_Training_RF_Locations = Unique_Training_RF_Locations;
    LHeiserAnalysis.retinalSources               = retinalSources;
    LHeiserAnalysis.retinal_preference           = retinal_preference;
    LHeiserAnalysis.SI                           = SI;
    
end