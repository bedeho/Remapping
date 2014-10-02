%
%  escape.m
%  Remapping
%
%  Created by Bedeho Mender on 23/09/14.
%  Copyright 2014 OFTNAI. All rights reserved.
%

% Platform specific escaping of uri
function escaped_str = escape(str)
    
    % windows , replace '\' with '\\'
    escaped_str = strrep(str, '\', '\\');
    
    %*nix, nothing needed
    % escaped_str = str;
end