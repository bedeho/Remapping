
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
        
        stimuli  = load([STIMULI_FOLDER stimulinames{i} filesep 'stim.mat']);
        type = stimuli.stimulitype;
        
        activityFile = [netDir filesep 'activity-' stimulinames{i} '.mat'];
        
        if strcmp(type,'StimuliControlTask'),
            receptivefield = AnalyzeStimuliControlTask(activityFile, stimulinames{i});
        elseif strcmp(type,'SaccadeControlTask'),
            receptivefield = AnalyzeSaccadeControlTask(activityFile, stimulinames{i});
        elseif strcmp(type,'KusonokiTesting'),
            receptivefield = AnalyzeStimuliControlTask(activityFile, stimulinames{i});
        else
            disp(['Unkonwn stimuli, skipping analysis...']);
        end
            
    end
    
    analysisSummary = 0;

end