
%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Analyze(netDir, stimulinames)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;
    
    % keep track of thesee
    stim_control_activity  = [];
    sacc_control_activity = [];
    stim_stimuli = [];
    sacc_stimuli = [];
    
    % Iterate stimuli
    for i = 1:length(stimulinames),
        
        % Load stimuli
        disp(['Loading stimuli: ' stimulinames{i}]);
        stimuliFile = [STIMULI_FOLDER stimulinames{i} filesep 'stim.mat'];
        stimuli  = load(stimuliFile);
        type = stimuli.stimulitype;
        
        % Load activity
        disp(['Loading activity: ' stimulinames{i}]);
        activityFile = [netDir filesep 'activity-' stimulinames{i} '.mat'];
        activity = load(activityFile);
        %activity = LoadActivity(activityFile);
        
        if strcmp(type,'StimuliControl'),
            
            disp('Doing stimuli control task analysis...');
            StimuliControl_Result = AnalyzeStimuliControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'StimuliControl_Result');
            
            stim_stimuli = stimuli;
            stim_control_activity = activity;
            stim_control_results = StimuliControl_Result;
            
        elseif strcmp(type,'SaccadeControl'),
            
            disp('Doing saccade control task analysis...');
            SaccadeControl_Result = AnalyzeSaccadeControlTask(activity, stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'SaccadeControl_Result');
            
            sacc_stimuli = stimuli;
            sacc_control_activity = activity;
            
        elseif strcmp(type,'DuhamelRemapping'),
            
            disp('Doing duhamel remapping task analysis...');
            DuhamelRemapping_Result = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemapping_Result');
            
        elseif strcmp(type,'DuhamelRemappingTrace'),
            
            disp('Doing duhamel remapping trace task analysis...');
            DuhamelRemappingTrace_Result = AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemappingTrace_Result');
            
        elseif strcmp(type,'DuhamelTruncation'),
            
            DuhamelTruncation_Result = AnalyzeDuhamelTruncation(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli);
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelTruncation_Result');
              
        elseif strcmp(type,'Kusonoki'),
            
            disp('Doing Kusonoki analysis...');
            [kusonokiSTIMAlignedAnalysis, kusonokiSACCAlignedAnalysis] = AnalyzeKusonoki(activity, stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'kusonokiSTIMAlignedAnalysis', 'kusonokiSACCAlignedAnalysis');
            
        elseif strcmp(type,'CLayerProbe'),
            
            disp('Doing C Layer Probe analysis...');
            [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity,  stimuli);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'CLabeProbe_Neurons_S', 'CLabeProbe_Neurons_V');
            
        else
            disp(['Unkonwn stimuli: ' num2str(stimulinames{i})]);
        end
            
    end
    
    % Plots
    [stmCtrlFigure, remScatFig, remTraceScatFig, kusonokiSACCFigure, kusonokiSTIMFigure, CLayerProbeFigure] = ThesisSimulationPlot(experimentFolder);
    
    % Stimuli Control
    saveas(stmCtrlFigure,[netDir filesep 'StimuliControl-summary.png']);
    saveas(stmCtrlFigure,[netDir filesep 'StimuliControl.eps']);
    close(stmCtrlFigure);

    %saveas(f,[netDir filesep 'DuhamelTruncation-summary.png']);
    %
    
    % Duhamel remapping trace
    saveas(remTraceScatFig,[netDir filesep name 'DuhamelRemappingTrace-summary.png']);
    saveas(remTraceScatFig,[netDir filesep name 'DuhamelRemappingTrace.eps']);
    close(remTraceScatFig);
    
    % Duhamel remaping
    saveas(remScatFig,[netDir filesep name 'DuhamelRemapping-summary.png']);
    saveas(remScatFig,[netDir filesep name 'DuhamelRemapping.eps']);
    close(remScatFig);
    
    %{
    
    % Kusonoki
    saveas(kusonokiSACCFigure,[netDir filesep 'Kusonoki-sacc-summary.png']);
    saveas(kusonokiSACCFigure,[netDir filesep 'Kusonoki-sacc.eps']);
    close(kusonokiSACCFigure);
    
    saveas(kusonokiSTIMFigure,[netDir filesep 'Kusonoki-summary.png']);
    saveas(kusonokiSTIMFigure,[netDir filesep 'Kusonoki.eps']);
    close(kusonokiSTIMFigure);
    
    % ClayerProbe
    
    saveas(CLayerProbeFigure,[netDir filesep 'CLayerProbe-summary.png']);
    saveas(CLayerProbeFigure,[netDir filesep 'CLayerProbe.eps']);
    close(CLayerProbeFigure);
    
    %}
end