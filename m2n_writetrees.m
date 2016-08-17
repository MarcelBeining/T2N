function tree = m2n_writetrees(params,tree,options,savepath)
% transforms the tree file into hoc code and also saves a interface file
% for correct node-section assignment
%
% options
% -cl: Cluster mode. Files are written to Server
% removed -d: Debug mode. Some measures are shown

if nargin < 3
    options = '';
end

sflag = 0;
if isstruct(tree)
    tree = {tree};
    sflag = 1;
end

if ~isfield(params,'path')
    params.path = regexprep(pwd,'\\','/');
else
    params.path = regexprep(params.path,'\\','/');
end
if isfield(params,'morphfolder')
    morphfolder = fullfile(params.path,params.morphfolder);
else
    errordlg('Please give morphfolder in params.morphfolder');
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

artflag = false(numel(tree),1);
for t=1:numel(tree)     % make neuron templates from trees and save/get minterface file
    artflag(t) = false;
    if ~isfield(tree{t},'artificial')
        [tree{t}, order] = sort_tree(tree{t},'-LO');
        if ~all(order==sort(order))
            orderchanged = true;
        end
    else
        artflag(t) = true;
    end
    nonameflag = false;
    if ~artflag(t) && isfield(params,'tname') && ischar(params.tname) && ~isempty(params.tname)
        treename = params.tname;
        countupflag = true;
    elseif artflag(t) && ~isfield(tree{t},'name')
        treename = tree{t}.artificial;
        countupflag = false;%true; % caution, maybe countup in artificial is necessary for correct assignment of records etc...check
        nonameflag = true;
    elseif isfield(tree{t},'name') % this is also true if name exist and it is artificial. good becaus
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
    if any(strfind(treename,'.'))
        badchars = badchars +numel(strfind(treename,'.'));
        treename(strfind(treename,'.')) = '_';
    end
    if (numel(treename) < 5 || ~strcmp(treename(1:5),'cell_'))
        treename = strcat('cell_',treename);
    end
    if countupflag
        treename = sprintf('%s_%d%d',treename,floor(t/10),rem(t,10));
    end
%     if strfind(options,'-cl')
%         [params.server.connect, answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/%s.hoc',nrn_morphfolder,treename));
%         fchk =  ~isempty(answer{1});
%     else
%         fchk = exist(fullfile(morphfolder,sprintf('%s.hoc',treename)),'file');
%     end
    
    oname = treename;
    neuron_template_tree (tree{t}, fullfile(morphfolder,sprintf('%s.hoc',treename)), [], '-m');
    if strfind(options,'-cl')   %transfer files to server
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s.hoc',oname)),sprintf('%s/%s.hoc',nrn_morphfolder,oname));
        pause(0.1)
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.dat',oname)),sprintf('%s/%s_minterf.dat',nrn_morphfolder,oname));
        pause(0.1)
        params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.mat',oname)),sprintf('%s/%s_minterf.mat',nrn_morphfolder,oname));
    end
    
    tree{t}.NID = treename;
    if nonameflag
        tree{t}.name = treename;
    end
end

if sflag
    tree = tree{1};
end
if ~all(artflag)
    if nargin < 4
        display('Please resave trees to have the NEURON ID in each tree')
        save_tree(tree);
    else
        save_tree(tree,savepath);
    end
end
if badchars > 0
    %     warndlg(sprintf('Caution! %d bad chars had to be removed or replaced from the tree names since they cause writing errors! Please be sure to not use "%%" and "." in the names',badchars),'Bad characters removed');
end
% if strfind(options,'-d')
%     tim = toc(tim);
%     fprintf(sprintf('Tree writing/reading time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
% end

if orderchanged && nargout == 0
    warndlg('Caution, the node order of some trees had to be changed! Sort your trees with "sort_tree" to obtain the correct results','Node order change!')
end
