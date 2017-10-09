% Run this code when you extracted the T2N zip file including this script at its final loation.
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

p = fileparts(mfilename('fullpath')); % get path to this script
addpath(genpath(p)); % add folder including subfolders to the Matlab path
fprintf('Added %s and its subdirectories to the Matlab path!\n',p);
