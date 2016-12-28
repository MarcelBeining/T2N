function t2n_init_modelfolders(folder)

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

mkdir(folder,'lib_mech')
mkdir(folder,'lib_custom')
mkdir(folder,'morphos')


% check for standard hoc files in the model folder and copy them if not existing
t2npath = which('t2n.m');  % get folder of t2n for copying files from it
if ~exist(fullfile(folder,'lib_genroutines'),'file')
    mkdir(folder,'lib_genroutines')
    display('non-existent folder lib_genroutines created')
end
if ~exist(fullfile(folder,'lib_genroutines/fixnseg.hoc'),'file')
    copyfile(fullfile(t2npath(1:end-6),'fixnseg.hoc'),fullfile(folder,'lib_genroutines/fixnseg.hoc'))
    display('fixnseg.hoc copied to model folder')
end
if ~exist(fullfile(folder,'lib_genroutines/genroutines.hoc'),'file')
    copyfile(fullfile(t2npath(1:end-6),'genroutines.hoc'),fullfile(folder,'lib_genroutines/genroutines.hoc'))
    display('genroutines.hoc copied to model folder')
end
if ~exist(fullfile(folder,'lib_genroutines/pasroutines.hoc'),'file')
    copyfile(fullfile(t2npath(1:end-6),'pasroutines.hoc'),fullfile(folder,'lib_genroutines/pasroutines.hoc'))
    display('pasroutines.hoc copied to model folder')
end

display('Please put your morphologies (.swc, .mtr, .neu etc.) into the "morphos" folder')
display('Please put mod files into the "lib_mech" folder')
display('Please put custom code (if needed) into the "lib_custom" folder and define with "neuron.custom" (see Documentation)')

