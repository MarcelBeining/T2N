function rename_nrnmech(newname,path)
% renames the nrnmech dll file  to the name specified in newname. Also
% deletes all .o and .c files created during dll compilation by NEURON
% INPUT
% newname: string with new name for dll file
% path (optional): path too the main folder of the model (containing the
% folder lib_mech)
%
% this function is part of the T2N package
% Copyright by Marcel Beining <marcel.beining@gmail.com>

if nargin < 1
    newname = 'nrnmech.dll';
end
if nargin < 2
    path = pwd;
end
origpth = pwd;

if isnumeric(newname)
    nam = sprintf('nrnmech_win%d.dll',newname);
elseif strcmp(newname(end-3:end),'.dll')
    nam = newname;
else
    nam = sprintf('%s.dll',newname);
end

if isempty(strfind(path,'lib_mech'))
    if ~exist(fullfile(path,'lib_mech'),'file')
        errordlg('No folder lib_mech exists')
        return
    else
       cd(fullfile(path,'lib_mech'))
    end
end
if ~strcmp('nrnmech.dll',nam)
    delete(nam)
    movefile('nrnmech.dll',nam)
end
fils = dir();
fils = fils(3:end);
fils = fils(cellfun(@(x) strcmp(x(end-1:end),'.o')|strcmp(x(end-1:end),'.c'),{fils.name}));
delete(fils.name)
cd(origpth)
display(sprintf('nrnmech.dll successfully renamed in %s ... and .o/.c files deleted',nam))