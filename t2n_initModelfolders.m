function t2n_initModelfolders(folder)
% This function initializes a new folder in which all necessary files and
% folders are created to create a compartmental model featuring T2N
%
% INPUTS
% folder                 (optional) path to the folder which should be
%                        initialized/created.
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 1
    folder = uigetdir(pwd,'Please give a folder where the model structure should be initialized');
end
if ~ischar(folder)
    errordlg('Input was no string')
    return
end
if ~exist(folder,'file')
    mkdir(folder)
end
% cd(folder)
if ~ exist(fullfile(folder,'lib_mech'),'dir')
    mkdir(folder,'lib_mech')
end
if ~ exist(fullfile(folder,'lib_custom'),'dir')
    mkdir(folder,'lib_custom')
end
if ~ exist(fullfile(folder,'morphos'),'dir')
    mkdir(folder,'morphos')
end


% check for standard hoc files in the model folder and copy them if not existing
t2npath = fileparts(which('t2n.m'));  % get folder of t2n for copying files from it
addpath(genpath(fullfile(t2npath,'src')))  % add src files to Matlab search path
if ~exist(fullfile(folder,'lib_genroutines'),'file')
    mkdir(folder,'lib_genroutines')
    display('non-existent folder lib_genroutines created')
end
if ~exist(fullfile(folder,'lib_genroutines/fixnseg.hoc'),'file')
    copyfile(fullfile(t2npath,'src','fixnseg.hoc'),fullfile(folder,'lib_genroutines/fixnseg.hoc'))
    display('fixnseg.hoc copied to model folder')
end
if ~exist(fullfile(folder,'lib_genroutines/genroutines.hoc'),'file')
    copyfile(fullfile(t2npath,'src','genroutines.hoc'),fullfile(folder,'lib_genroutines/genroutines.hoc'))
    display('genroutines.hoc copied to model folder')
end
if ~exist(fullfile(folder,'lib_genroutines/pasroutines.hoc'),'file')
    copyfile(fullfile(t2npath,'src','pasroutines.hoc'),fullfile(folder,'lib_genroutines/pasroutines.hoc'))
    display('pasroutines.hoc copied to model folder')
end

display('Please put your morphologies (.swc, .mtr, .neu etc.) into the "morphos" folder')
display('Please put mod files into the "lib_mech" folder')
display('Please put custom code (if needed) into the "lib_custom" folder and define with "neuron.custom" (see Documentation)')

