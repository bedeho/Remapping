
%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function analysisSummary = Analyze(netDir, stimulinames)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;

    for i = 1:length(stimulinames),
        
        stimuliFile = [STIMULI_FOLDER stimulinames{i} filesep 'stim.mat'];
        stimuli  = load(stimuliFile);
        type = stimuli.stimulitype;
        dt = stimuli.dt;
        
        activityFile = [netDir filesep 'activity-' stimulinames{i} '.mat'];
        
        if strcmp(type,'StimuliControlTask'),
            
            [baselineResponse, stim_response, location, foundOnset, foundOffset, latencyTimeStep, durationTimeStep, neuronResponse] = AnalyzeStimuliControlTask(activityFile, stimuliFile);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , ...
                    'baselineResponse', ...
                    'stim_response', ...
                    'location', ...
                    'foundOnset', ...
                    'foundOffset', ...
                    'latencyTimeStep', ...
                    'durationTimeStep');
                            
            R_N = size(latencyTimeStep,2);                
            f = figure;
            imagesc(neuronResponse);
            ylabel('Neuron');
            xlabel('Time');
            hold on;
            plot(latencyTimeStep,1:R_N,'wo');
            saveas(f,[netDir filesep stimulinames{i} '.png']);
            
        elseif strcmp(type,'SaccadeControlTask'),
            
            saccade_response = AnalyzeSaccadeControlTask(activityFile, stimuliFile);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'saccade_response');
            
        elseif strcmp(type,'KusonokiTesting'),
            
            [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activityFile, stimulinames{i});
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'kusonokiSTIMAlignedAnalysis', 'kusonokiSACCAlignedAnalysis');
            
        else
            disp(['Unkonwn stimuli: ' num2str(stimulinames{i})]);
        end
            
    end
    
    analysisSummary = 0;
end