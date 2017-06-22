function tree = t2n_writeTrees(tree,params,savepath,options)
% This transforms the tree file into hoc code and also saves a interface file
% for correct node-section assignment
%
% INPUTS
% tree              tree cell array with morphologies (see documentation)
% params            t2n parameter structure (see documentation)
% savepath          (optional) if this file destination string is given, the
%                   function does not have to ask for it via gui
% options           string with optional arguments (can be concatenated):
%                   -d: Debug mode. The duration of writing hocs is shown
%                   -w: Waitbar showing the progress
%                   %deactived for the moment: -cl: Cluster mode. Files are written to Server
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 4
    options = '';
end


sflag = 0;
if isstruct(tree)
    tree = {tree};
    sflag = 1;
end
if size(tree,1) ~= numel(tree)  % check for correct 1st dimension 
    tree = tree';
end
if ~isfield(params,'path')
    params.path = regexprep(pwd,'\\','/');
else
    params.path = regexprep(params.path,'\\','/');
end
if isfield(params,'morphfolder')
    morphfolder = fullfile(params.path,params.morphfolder);
elseif exist(fullfile(params.path,'morphos'),'dir')
    morphfolder = fullfile(params.path,'morphos');
else
    errordlg('Please give morphology folder in "params.morphfolder" or create the standard folder "morphos" using "t2n_init_modelfolders"!');
    return
end

if strfind(options,'-cl')
    if isfield(params,'morphfolder')
        nrn_morphfolder = fullfile(params.server.clpath,params.morphfolder);
    else
        nrn_morphfolder = nrn_exchfolder;
    end
    
    sshfrommatlabissue(params.server.connect,sprintf('mkdir -p %s',nrn_morphfolder));
else
    nrn_morphfolder = morphfolder;
end
nrn_morphfolder = regexprep(nrn_morphfolder,'\\','/');

if ~exist(morphfolder,'file')
    mkdir(morphfolder);
end

orderchanged = false;
badchars = 0;

artflag = cellfun(@(x) isfield(x,'artificial'),tree);  % get boolean for trees being artificial
indartflag = find(artflag);   % get indices
[~,ia] = unique(cellfun(@(x) x.artificial,tree(artflag),'UniformOutput',0));  % find artificial cells that are of the same type
indWrite = cat(1,find(~artflag),indartflag(ia));  % only write trees that are not artificial and one artificial tree of each type
% tree = cat(1,tree(~artflag),tree(indartflag(ia)));  % only write trees that are not artificial and one artificial tree of each type

tim = tic;
if ~isempty(strfind(options,'-w'))
    w = waitbar(0,'Trees are transformed to hoc, please wait...');
end
for t=1:numel(tree)     % make neuron templates from trees and save/get minterface file
    if ~artflag(t)
        [tree{t}, order] = sort_tree(tree{t},'-LO');
        if ~all(order==sort(order))
            orderchanged = true;
        end
    end
    nonameflag = false;
    if ~artflag(t) && isfield(params,'tname') && ischar(params.tname) && ~isempty(params.tname)
        treename = params.tname;
        countupflag = true;
    elseif artflag(t) && ~isfield(tree{t},'name')
        treename = tree{t}.artificial;
        countupflag = false;
        nonameflag = true;
    elseif isfield(tree{t},'name') % this is also true if name exist and it is artificial
        treename = tree{t}.name;
        countupflag = false;
    else
        treename = 'Tree';
        countupflag = true;
        nonameflag = true;
    end
    if any(strfind(treename,'%'))
        badchars = badchars +numel(strfind(treename,'%'));
        treename(strfind(treename,'%')) = [];
    end
    if any(strfind(treename,'-'))
        badchars = badchars +numel(strfind(treename,'%'));
        treename(strfind(treename,'-')) = '_';
    end
    if any(strfind(treename,'.'))
        badchars = badchars +numel(strfind(treename,'.'));
        treename(strfind(treename,'.')) = '_';
    end
    if (numel(treename) < 5 || ~strcmp(treename(1:5),'cell_'))
        treename = strcat('cell_',treename);
    end
    if countupflag
        treename = sprintf('%s_%d',treename,t);
    end
%     if strfind(options,'-cl')
%         [params.server.connect, answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/%s.hoc',nrn_morphfolder,treename));
%         fchk =  ~isempty(answer{1});
%     else
%         fchk = exist(fullfile(morphfolder,sprintf('%s.hoc',treename)),'file');
%     end
if any(t == indWrite)  % inly rewrite artificial trees once
    oname = treename;
    neuron_template_tree (tree{t}, fullfile(morphfolder,sprintf('%s.hoc',treename)), '-m');
    
    if strfind(options,'-cl')   %transfer files to server
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s.hoc',oname)),sprintf('%s/%s.hoc',nrn_morphfolder,oname));
        pause(0.1)
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.dat',oname)),sprintf('%s/%s_minterf.dat',nrn_morphfolder,oname));
        pause(0.1)
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.mat',oname)),sprintf('%s/%s_minterf.mat',nrn_morphfolder,oname));
    end
end
    tree{t}.NID = treename;
    if nonameflag || badchars > 0
        tree{t}.name = treename;
    end
    if ~isempty(strfind(options,'-w'))
        waitbar(t/numel(tree),w);
    end
end

if sflag
    tree = tree{1};
end
if ~all(artflag)
    if nargin < 3
        display('Please resave trees to have the NEURON ID in each tree')
        save_tree(tree);
    else
        save_tree(tree,savepath);
    end
end
if ~isempty(strfind(options,'-w'))
    close(w)
end
if badchars > 0
        warndlg(sprintf('Caution! %d bad chars had to be removed or replaced from the tree names since they cause writing errors! Please be sure to not use "%%" and "." in the names',badchars),'Bad characters removed');
end
if strfind(options,'-d')
    tim = toc(tim);
    fprintf(sprintf('Tree hoc writing time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
end

if orderchanged && nargout == 0
    warndlg('Caution, the node order of some trees had to be changed! Sort your trees with "sort_tree" to obtain the correct results','Node order change!')
end
