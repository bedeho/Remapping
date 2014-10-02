%
%  declareGlobalVars.m
%  Remapping
%
%  Created by Bedeho Mender on 19/05/13.
%  Copyright 2013 OFTNAI. All rights reserved.
%

function declareGlobalVars()

	global base;
	base = ['C:' filesep 'Users' filesep 'bedeho' filesep 'Documents' filesep 'GitHub' filesep 'Remapping' filesep]; % must have trailing slash

	global EXPERIMENTS_FOLDER;
	EXPERIMENTS_FOLDER = [base 'Experiments' filesep];
    
    global STIMULI_FOLDER;
	STIMULI_FOLDER = [base 'Stimuli' filesep];
    
    global THESIS_FIGURE_PATH;
    THESIS_FIGURE_PATH = ['C:' filesep 'Users' filesep 'bedeho' filesep 'Documents' filesep 'GitHub' filesep 'Remapping' filesep 'Thesisfigures' filesep ''];
    
    