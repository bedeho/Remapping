
%
%  Analyze.m
%  Analyze
%
%  Created by Bedeho Mender on 11/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function Sprattling_Analyze(experimentFolder, netDir, stimulinames, network)

    % Import global variables
    declareGlobalVars();
    global STIMULI_FOLDER;
    
    % keep track of thesee
    stim_control_activity  = [];
    sacc_control_activity = [];
    stim_stimuli = [];
    sacc_stimuli = [];
    
    isUntrained = strcmp(network,'BlankNetwork');
    
    % Iterate stimuli
    for i = 1:length(stimulinames),
        
        % Load stimuli
        disp(['Loading stimuli: ' stimulinames{i}]);
        
        % OLD GLOBAL
        %stimuliFile = [STIMULI_FOLDER stimulinames{i} filesep 'stim.mat'];
        
        %LOCAL
        stimuliFile = [experimentFolder 'STIM-' stimulinames{i} filesep 'stim.mat'];
        
        stimuli  = load(stimuliFile);
        type = stimuli.stimulitype;
        
        % Load activity
        disp(['Loading activity: ' stimulinames{i}]);
        activityFile = [netDir filesep 'activity-' stimulinames{i} '.mat'];
        activity = load(activityFile);
        %activity = LoadActivity(activityFile);
        
        if strcmp(type,'StimuliControl'),
            
            disp('Doing stimuli control task analysis...');
            [StimuliControl_Result, Decoded_ReceptiveFieldsLocations] = Sprattling_AnalyzeStimuliControlTask(activity, stimuli);
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
            
       elseif strcmp(type,'CLayerProbe'),
            
            disp('Doing C Layer Probe analysis...');
            [CLabeProbe_Neurons_S, CLabeProbe_Neurons_V] = AnalyzeCLayerProbe(activity,  stimuli);
            
            % We also save limits of input to C map for convenicence
            % uring plotting
            R_max = stimuli.R_eccentricity;
            S_max = stimuli.S_eccentricity;
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'CLabeProbe_Neurons_S', 'CLabeProbe_Neurons_V', 'R_max', 'S_max');
            
        elseif strcmp(type,'DuhamelRemappingTrace') && ~isUntrained,
            
            disp('Doing duhamel remapping trace task analysis...');
            
            [DuhamelRemappingTrace_Result] = Sprattling_AnalyzeDuhamelRemapping(activity, stimuli, stim_control_activity.R_firing_history, stim_stimuli, sacc_control_activity.R_firing_history, sacc_stimuli, Decoded_ReceptiveFieldsLocations);
            
            save([netDir filesep 'analysis-' stimulinames{i} '.mat'] , 'DuhamelRemappingTrace_Result');
            
        else
            disp(['Unsupported stimuli: ' num2str(stimulinames{i})]);
        end
            
    end
    
   %% Plots
   
   % Stimuli Control
   if(~isUntrained),
        f = figure;   
        hist(Decoded_ReceptiveFieldsLocations, 40);
        title('Receptive Field Location Distribution');
        saveas(f,[netDir filesep 'StimuliControl-summary.png']);
        saveas(f,[netDir filesep 'StimuliControl.eps'], 'epsc');
        close(f);
   end
        
    
    %{
    % Stimuli Control
    if(remTraceScatFig),
        saveas(stmCtrlFigure,[netDir filesep 'StimuliControl-summary.png']);
        saveas(stmCtrlFigure,[netDir filesep 'StimuliControl.eps'], 'epsc');
        close(stmCtrlFigure);
    end
    
    % Duhamel remapping trace
    if(remTraceScatFig),
        saveas(remTraceScatFig,[netDir filesep 'DuhamelRemappingTrace-summary.png']);
        saveas(remTraceScatFig,[netDir filesep 'DuhamelRemappingTrace.eps'], 'epsc');
        close(remTraceScatFig);
    end
    
    % Duhamel remaping
    if(remScatFig),
        saveas(remScatFig,[netDir filesep 'DuhamelRemapping-summary.png']);
        saveas(remScatFig,[netDir filesep 'DuhamelRemapping.eps'], 'epsc');
        close(remScatFig);
    end
    
    % LHeiser
    if(LHeiserFig),
        saveas(LHeiserFig,[netDir filesep 'LHeiser-summary.png']);
        saveas(LHeiserFig,[netDir filesep 'LHeiser.eps'], 'epsc');
        close(LHeiserFig);
    end
    
    % Kusonoki
    if(kusonokiSACCFig),
        saveas(kusonokiSACCFig,[netDir filesep 'Kusonoki-sacc-summary.png']);
        saveas(kusonokiSACCFig,[netDir filesep 'Kusonoki-sacc.eps'], 'epsc');
        close(kusonokiSACCFig);
    end
    
    if(kusonokiSTIMFig),
        saveas(kusonokiSTIMFig,[netDir filesep 'Kusonoki-summary.png']);
        saveas(kusonokiSTIMFig,[netDir filesep 'Kusonoki.eps'], 'epsc');
        close(kusonokiSTIMFig);
    end
    
    % ClayerProbe
    if(CLayerProbeFigure),
        saveas(CLayerProbeFigure,[netDir filesep 'CLayerProbe-summary.png']);
        saveas(CLayerProbeFigure,[netDir filesep 'CLayerProbe.eps']);
        close(CLayerProbeFigure);
    end
    %}
end