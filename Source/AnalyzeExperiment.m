
%
%  AnalyzeExperiment.m
%  Remapping
%
%  Created by Bedeho Mender on 24/06/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function AnalyzeExperiment()

            
                                % Move each network to new folder
                                listing = dir(experimentFolder); 
                                for d = 1:length(listing),

                                    % We looking for networks files
                                    subsim_name = listing(d).name;

                                    if ~listing(d).isdir && ~isempty(findstr(subsim_name,'Network')),

                                        % Make dir name
                                        subsim_dir = [simulationFolder filesep subsim_name];

                                        % Make dir
                                        mkdir();

                                        % Move file into dir
                                        movefile([simulationFolder filesep subsim_name],\);
                                    end
                                end


                                % Testing & Analysis: Later iterate all subfolder that may have been created 
                                Analyze(simulationFolder, testingStimuliFile);
    
end