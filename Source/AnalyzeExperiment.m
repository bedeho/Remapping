
%
%  AnalyzeExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 24/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function AnalyzeExperiment(experiment, stimulinames, trainingStimuli)

    % Experiment name 
    if nargin < 1,
        experiment = 'sprattling_visual_learning_bigepoch20';
    end
    
    % Stimuli names
    if(nargin < 3)
        
        trainingStimuli = 'basic-Training';
        
        %trainingStimuli = 'basic-Training_LHeiser';
        
        if(nargin < 2)
            
            %{
            stimulinames = {'basic-StimuliControl', ...
                            'basic-SaccadeControl', ...
                            'basic-CLayerProbe', ...
                            'basic-DuhamelRemapping', ...
                            'basic-DuhamelRemappingTrace', ... 
                            'basic-DuhamelTruncation', ...
                            'basic-Kusonoki'};
            %}
            
            %, 'basic-CLayerProbe',
            
            stimulinames = {'basic-StimuliControl', 'basic-SaccadeControl', 'basic-DuhamelRemappingTrace'};
            
        end
    end
    
    % Import global variables
    declareGlobalVars();
    
    global STIMULI_FOLDER;
    global EXPERIMENTS_FOLDER;
    global base;

    experimentFolder = [EXPERIMENTS_FOLDER experiment filesep];
    
    % Iterate simulations in this experiment folder
    listing = dir(experimentFolder); 
        
    % Find an example of simulation directory to extract column names- HORRIBLY CODED
    for d = 1:length(listing),

        simulation = listing(d).name;

        if listing(d).isdir && ~any(strcmp(simulation, {'Filtered', 'Images', '.', '..'})) && ~strcmp(simulation(1:5), 'STIM-'),
            [parameters, nrOfParams] = getParameters(simulation);
            break;
        end
    end

    % Save results for summary
    filename = [experimentFolder 'Summary.html']; % Make name that is always original so we dont overwrite old summary which is from previous xGridCleanup run of partial results from this same parameter search
    fileID = fopen(filename, 'w'); % did note use datestr(now) since it has string
    
    % HTML header with javascript/css for JQUERY table plugin
    fprintf(fileID, '<html><head>\n');
    fprintf(fileID, '<style type="text/css" title="currentStyle">\n');
                             
    fprintf(fileID, ['@import "' escape([base 'Source' filesep 'DataTables-1.8.2' filesep 'media' filesep 'css' filesep 'demo_page.css"']) ';\n']);
	fprintf(fileID, ['@import "'  escape([base 'Source' filesep 'DataTables-1.8.2' filesep 'media' filesep 'css' filesep 'demo_table.css"']) ';\n']);
	fprintf(fileID, '</style>\n');
	fprintf(fileID, ['<script type="text/javascript" language="javascript" src="' escape([base 'Source' filesep 'DataTables-1.8.2' filesep 'media' filesep 'js' filesep 'jquery.js']) '"></script>\n']);
	fprintf(fileID, ['<script type="text/javascript" language="javascript" src="' escape([base 'Source' filesep 'DataTables-1.8.2' filesep 'media' filesep 'js' filesep 'jquery.dataTables.js']) '"></script>\n']);
	fprintf(fileID, '<script type="text/javascript" charset="utf-8">\n');
	fprintf(fileID, '$(document).ready(function() { $("#example").dataTable();});\n');
	fprintf(fileID, '</script>\n');
    fprintf(fileID, '</head>\n');
    fprintf(fileID, '<body>\n');
    
    % HTML Title
    fprintf(fileID, '<h1>%s - %s</h1>\n', experiment, datestr(now));
    
    % HTML table
    fprintf(fileID, '<table id="example" class="display" cellpadding="10" style="border: solid 1px">\n');
    fprintf(fileID, '<thead><tr>');
    fprintf(fileID, '<th>Simulation</th>');
    fprintf(fileID, '<th>Network</th>');
    
    for i = 1:nrOfParams,
        fprintf(fileID, ['<th>' parameters{i,1} '</th>']);
    end
    
    fprintf(fileID, '<th>Summary</th>');
    
    stimuli_types = cell(1,length(stimulinames))';
    for i = 1:length(stimulinames),
        
        fprintf(fileID, ['<th>' stimulinames{i} '</th>']);
        stimuli             = load([STIMULI_FOLDER stimulinames{i} filesep 'stim.mat']);
        stimuli_types{i}    = stimuli.stimulitype;
        stimuli_file{i}     = [experimentFolder 'STIM-' stimulinames{i}];
        
    end
    
    trainingStimuliFile = [experimentFolder 'STIM-' trainingStimuli];
    
    counter = 1;
    
    format('short');
    
    fprintf(fileID, '<tbody>\n');
    for d = 1:length(listing),

        % We are only looking for directories, but not the
        % 'Filtered' directory, since it has filtered output
        simulation = listing(d).name;

        if listing(d).isdir && ~any(strcmp(simulation, {'Filtered', 'Images', '.', '..'})) && ~strcmp(simulation(1:5), 'STIM-'),
            
            disp(['******** Simulation ' num2str(counter) ' ********']); % ' out of ' num2str((nnz([listing(:).isdir]) - 2)) '
            counter = counter + 1;
            
            % Iterate simulations in this experiment folder
            simulationFolder = [experimentFolder simulation];
            simulationListing = dir(simulationFolder); 

            for s=1:length(simulationListing),
                
                network = simulationListing(s).name; 
                
                if simulationListing(s).isdir && ~any(strcmp(network, {'Training', '.', '..'})),
                    
                    % color all these rows in same color
                    
                    netDir = [simulationFolder filesep network];
                    
                    %netDirRelative = [simulation filesep network];
                    %trDir = [experimentFolder simulation filesep 'Training'];
                    
                    % Do analysis
                    %Analyze(experimentFolder, netDir, stimulinames);
                    Sprattling_Analyze(experimentFolder, netDir, stimulinames, network);
                    
                    % Start row
                    fprintf(fileID, '<tr>');

                    % Name
                    fprintf(fileID, '<td> %s </td>\n', simulation);

                    % Network
                    fprintf(fileID, '<td> %s </td>\n', network);

                    % Parameters
                    parameters = getParameters(simulation);

                    for i = 1:nrOfParams,
                        fprintf(fileID, ['<td> ' parameters{i,2} ' </td>\n']);
                    end
                    
                    % Summary
                    fprintf(fileID, '<td>');
                    %fprintf(fileID, '<img src="%s" width="350px" height="350px"/>\n', [netDir filesep 'summary.png']);
                    outputButton('Training', ['matlab:viewNeuronDynamics(\\''' escape(escape([simulationFolder filesep 'activity-' trainingStimuli '.mat'])) '\\'',\\''' escape(escape(trainingStimuliFile)) '\\'',\\''' escape(escape([netDir filesep network '.mat'])) '\\'',\\''' escape(escape([netDir filesep 'analysis-basic-CLayerProbe.mat'])) '\\'')']);
                    fprintf(fileID, '</td>');
                    
                    % Stimuli
                    for i = 1:length(stimulinames),
                        
                        fprintf(fileID, '<td>');
                        
                        % Image - only if it exists!
                        imageFile = [netDir filesep stimuli_types{i} '-summary.png'];
                        
                        if exist(imageFile, 'file'),                        
                            fprintf(fileID, '<img src="%s" width="250px" height="250px"/></br>\n', imageFile); %400
                        end
                        
                        % Button
                        outputButton('Activity', ['matlab:viewNeuronDynamics(\\''' escape(escape([netDir filesep 'activity-' stimulinames{i} '.mat'])) '\\'',\\''' escape(escape(stimuli_file{i})) '\\'',\\''' escape(escape([netDir filesep network '.mat'])) '\\'',\\''' escape(escape([netDir filesep 'analysis-basic-CLayerProbe.mat'])) '\\'')']);
                        
                        fprintf(fileID, '</td>');
                    end

                    fprintf(fileID, '\n\n');
                
                end
            end
            
        end
    end

    fprintf(fileID, '</table></body></html>');
    fclose(fileID);
     
    % Thanks to DanTheMans excellent advice, now we dont
    % have to fire up terminal insessantly
    
    disp('Analysis completed...');
    system('stty echo');
    
    web(filename);
    
    function outputButton(title, action)
        fprintf(fileID, ['<input type="button" value="' title '" onclick="document.location=''' action '''"/></br>\n']);
    end

    function [parameters, nrOfParams] = getParameters(sim)

        % is there more than one paramter, handle it in special case if so.
        if ~isempty(strfind(sim, '-')),
            columns = strsplit(sim, '-');
            nrOfParams = length(columns);
            
            parameters = cell(nrOfParams,2);

            for p = 1:nrOfParams,
                pair = strsplit(char(columns(p)),'='); % columns(p) is a 1x1 cell
                parameters{p,1} = char(pair(1));
                parameters{p,2} = char(pair(2));
            end
        
        else
            
            if ~isempty(strfind(sim, '=')),
                
                nrOfParams = 1;
                
                pair = strsplit(sim,'=');
                parameters = cell(nrOfParams,2);
                parameters{1,1} = char(pair(1));
                parameters{1,2} = char(pair(2));
            else
                parameters = [];
                nrOfParams = 0;
            end
            
        end
    end
end