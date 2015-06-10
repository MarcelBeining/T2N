function tree = m2n_writetrees(params,tree,options)

if nargin < 3
    options = '';
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
if strfind(options,'-d')
    tim = tic;
end

for t=1:numel(tree)     % make neuron templates from trees and save/get minterface file
    artflag = false;
    if ~isfield(tree{t},'artificial')
        [tree{t}, order] = sort_tree(tree{t},'-LO');
        if ~all(order==sort(order))
            orderchanged = true;
        end
    else
        artflag = true;
    end
    if isfield(params,'tname') && ischar(params.tname) && ~isempty(params.tname)
        treename = params.tname;
        tflag = true;
    elseif artflag && ~isfield(tree{t},'name')
        treename = tree{t}.artificial;
        tflag = true;
    elseif isfield(tree{t},'name')
        treename = tree{t}.name;
        tflag = false;
    else
        treename = 'Tree';
        tflag = true;
    end
    if any(strfind(treename,'%'))
        badchars = badchars +numel(strfind(treename,'%'));
        treename(strfind(treename,'%')) = [];
    end
    if any(strfind(treename,'.'))
        badchars = badchars +numel(strfind(treename,'.'));
        treename(strfind(treename,'.')) = '_';
    end
    
    treename = strcat('cell_',treename);
    if tflag
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
end

display('Please resave trees to have the NEURON ID in each tree')
save_tree(tree);


if badchars > 0
    %     warndlg(sprintf('Caution! %d bad chars had to be removed or replaced from the tree names since they cause writing errors! Please be sure to not use "%%" and "." in the names',badchars),'Bad characters removed');
end
if strfind(options,'-d')
    tim = toc(tim);
    fprintf(sprintf('Tree writing/reading time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
end

if orderchanged
    warndlg('Caution, the node order of some trees had to be changed! Sort your trees with "sort_tree" to obtain the correct results','Node order change!')
end