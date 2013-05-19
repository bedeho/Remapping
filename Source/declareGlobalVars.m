%
%  declareGlobalVars.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function declareGlobalVars()

	global base;
	base = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Projects/Remapping/'; % must have trailing slash

	global EXPERIMENTS_FOLDER;
	EXPERIMENTS_FOLDER = [base 'Experiments' filesep];
    
    global THESIS_FIGURE_PATH;
    THESIS_FIGURE_PATH = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/mender/Dphil/Thesis/figures/';
    
    