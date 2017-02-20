function [out, origminterf,params,tree] = t2n(tree,params,neuron,options)
% function t2n ("Trees toolbox to Neuron") to generate and execute a NEURON
% simulation with the morphologies in tree, and parameters in the structure 
% params and neuron
% The output-file(s) of the NEURON function are read and transferred
% into the output variable out
%
% options:
%   -w waitbar
%   -d Debug mode (NEURON is opened and some parameters are set)
%   -q quiet mode -> suppress output and suppress NEURON instance(s) opening
%   -m let T2N recompile the nrnmech.dll. Useful if a mod file was modified. 
%      For safety of compiled dlls, this option does not work when an explicit 
%      name of a dll was given via params.nrnmech!
%   -cl cluster mode -> files are prepared to be executed on a cluster.
%
% This code originated from an idea of Johannes Kasper, a former group-member
% in the Lab of Hermann Cuntz, Frankfurt.
%
% Copyright by marcel.beining@gmail.com, June 2015

%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interf_file = 'neuron_runthis.hoc'; % the name of the main hoc file which will be written

%% check options and paths

t2npath = fileparts(which('t2n.m'));

if ~isfield(params,'path')
    params.path = regexprep(pwd,'\\','/');
    display('No standard path was given. Current folder is used instead');
else
    params.path = regexprep(params.path,'\\','/');
end

if nargin < 4 || isempty(options)
    options = '';
end
if ~isempty(strfind(options,'-d'))
    debug = 1;
else
    debug = 0;
end
if strcmpi(params.path(end),'\')
    params.path = params.path(1:end-1);
end

if ~isempty(strfind(options,'-cl')) % server mode (not fully implemented)
    nrn_path = params.server.modelfolder;
    if ~isfield(params.server,'walltime') || numel(params.server.walltime) ~=3
        params.server.walltime = [5 0 0];
        warning('Server walltime not specified correctly (1 by 3 vector in params.server.walltime). Walltime set to 5 hours');
    end
    if ~isfield(params.server,'memory') || numel(params.server.memory) ~=1
        params.server.memory = 1;
        warning('Max memory per node not specified correctly (1 scalar [GB] in params.server.memory). Max memory set to 1GB');
    end
    params.server.softwalltime = sum(params.server.walltime .* [3600,60,1])-30; %subtract 30 sec
    params.server.softwalltime = floor([params.server.softwalltime/3600,rem(params.server.softwalltime,3600)/60,rem(params.server.softwalltime,60)]);
    params.server.envstr = '';
    params.server.qfold = '';
    if isfield(params.server,'env') && isstruct(params.server.env)
        envnames = fieldnames(params.server.env);
        for v = 1:numel(envnames)
            params.server.envstr = [params.server.envstr,'setenv ',envnames{v},' ',params.server.env.(envnames{v}),'; '];
        end
        if isfield(params.server.env,'SGE_ROOT') && isfield(params.server.env,'SGE_ARCH')
            params.server.qfold = sprintf('%s/bin/%s/',params.server.env.SGE_ROOT,params.server.env.SGE_ARCH);
        end
    end
else
    nrn_path = params.path;
end
if strcmp(nrn_path(end),'/') % remove "/" from path
    nrn_path = nrn_path(1:end-1);
end


%% initialize basic variables
noutfiles = 0;          % counter for number of output files
readfiles = cell(0);    % cell array that stores information about output files

%% check other input

if nargin < 2 || isempty(params)
    if debug == 1
        %% individual parameters structure for debug
        params.openNeuron = true;
        params.nseg = 'dlambda';
        params.tstop = 200;
        params.dt = 0.025;
        params.accuracy = 0;
        params.prerun = false;
    else
        params = [];
    end
end

if nargin < 3 || isempty(neuron) % if no neuron structure was given, create standard parameter set
    if debug ~= 1
        warning('No input about what to do was given! Standard test (HH + rectangle somatic stimulation) is used')
    end
    neuron.pp{1}.IClamp = struct('node', 1, 'del',100,'dur',50,'amp',0.05);
    neuron.record{1}.cell = struct('node',1 ,'record', 'v');
    neuron.APCount{1} = {1,-30}; % {node, tresh}
    neuron.mech.all.pas = [];
    neuron.mech.soma.hh = [];
    neuron.mech.axon.hh = [];
end

outoptions.nocell = false;
if isstruct(neuron)         % transform structure neuron to cell neuron
    if numel(neuron) == 1
        neuron = {neuron};
        outoptions.nocell = true;
    else
        neuron = arrayfun(@(y) {y},neuron);
    end
end
out = cell(numel(neuron),1);

% check if trees are there and composed within a cell array
if nargin < 1 || isempty(tree)
    error('No tree specified in input')
end
if iscell(tree) && iscell(tree{1})
    tree = tree{1};
elseif isstruct(tree)
    tree = {tree};
end

% check for morphology folder
if isfield(params,'morphfolder')
    morphfolder = fullfile(params.path,params.morphfolder);
elseif exist(fullfile(nrn_path,'morphos'),'dir')
    morphfolder = fullfile(nrn_path,'morphos');
    params.morphfolder = 'morphos';
else
    error('Please give morphology folder in params.morphfolder or create standard morphology folder "morphos"');
%     return
end
morphfolder = regexprep(morphfolder,'\\','/');

if ~isfield(params,'neuronpath')
    params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';  % default path to neuron exe
end

if strfind(options,'-cl') % server mode, not fully implemented
    if ~isfield(params,'server')
        error('No access data provided for Cluster server. Please specify in params.server')
%         return
    else
        if isfield(params.server,'connect')
            
        else
            if  exist('sshj-master/sshj-0.0.0-no.git.jar','file')% exist('ganymed-ssh2-build250/ganymed-ssh2-build250.jar','file')
                sshjfolder = fileparts(which('sshj-master/sshj-0.0.0-no.git.jar'));
                javaaddpath(fullfile(sshjfolder,'sshj-0.0.0-no.git.jar'));%javaaddpath(which('ganymed-ssh2-build250/ganymed-ssh2-build250.jar'));
                javaaddpath(fullfile(sshjfolder,'bcprov-ext-jdk15on-156.jar'));
                javaaddpath(fullfile(sshjfolder,'slf4j-1.7.23'));
            else
                try
                    sshfrommatlabinstall(1)
                catch
                    error('Could not find the ganymed ssh zip file')
%                     return
                end
            end
            params.server.connect = sshfrommatlab(params.server.user,params.server.host,params.server.pw);
        end
        if ~isfield(params.server,'modelfolder')
            %            params.server.modelfolder = '~';
            %            warning('No Path on the Server specified. Root folder will be used')
            error('No folder on Server specified, please specify under params.server.modelfolder')
%             return
        end
    end
end


% check for exchange folder (folder where files between Matlab and NEURON
% are exchanged)
if isfield(params,'exchfolder')
    exchfolder = fullfile(params.path,params.exchfolder);
    if strfind(options,'-cl')
        nrn_exchfolder = fullfile(params.server.modelfolder,params.exchfolder);
    else
        nrn_exchfolder = exchfolder;
    end
else
    exchfolder = fullfile(params.path,'t2n_exchange');
    if strfind(options,'-cl')
        nrn_exchfolder = fullfile(params.server.modelfolder,'t2n_exchange');
    else
        nrn_exchfolder = exchfolder;
    end
    params.exchfolder = 't2n_exchange';
end
nrn_exchfolder = regexprep(nrn_exchfolder,'\\','/');


% check for several standard parameter and initialize default value if not set
if ~isfield(params,'cvode')
    params.cvode = false;
end
if ~isfield(params,'use_local_dt')
    params.use_local_dt = 0;
end
if ~isfield(params,'openNeuron') || ~isempty(strfind(options,'-q'))
    params.openNeuron = false;
end
if ~isfield(params,'nseg') || strcmpi(params.nseg,'d_lambda')
    params.nseg = 'dlambda';
    display('Number of segments or nseg rule not set in params.nseg. Dlambda rule will be applied')
end
if ~isfield(params,'tstart')
    params.tstart = 0;
end
if ~isfield(params,'tstop')
    params.tstop = 200;
    display('Tstop not defined in params.tstop. Default value of 200 ms is applied.')
end
if ~isfield(params,'dt')
    params.dt = 0.025;
    if ~params.cvode
        display('Time step not defined in params.dt. Default value of 0.025 ms is applied.')
    end
end
if ~isfield(params,'accuracy')
    params.accuracy = 0;
end
if ~isfield(params,'skiprun')
    params.skiprun = false;
end
if ~isfield(params,'q10')
    params.q10 = false;
end
if ~isfield(params,'prerun')
    params.prerun = false;
end


% Check if NEURON software exists at the given path
if ~isempty(strfind(options,'-cl'))
    [~, outp] = sshfrommatlabissue(params.server.connect,'module avail');
    params.server.neuron = regexpi(outp.StdErr,'neuron/\d{1,2}\.\d{1,2}\s','match');  % available modules are reported to errorStream..dunno why
%     params.server.envstr = [params.server.envstr, sprintf('module load %s; ',outp{1})];  % load first found neuron version
    display(params.server.neuron{1})
elseif exist(params.neuronpath,'file') ~= 2
    if isempty(strfind(params.neuronpath,'.exe') && exist(fullfile(params.neuronpath,'nrniv.exe'),'file') == 2) % path only points to folder, not to exe
        params.neuronpath = fullfile(params.neuronpath,'nrniv.exe');
    else
        error('No NEURON software (nrniv.exe) found under "%s"\nPlease give correct path using params.neuronpath',params.neuronpath);
%         return
    end
end

% check for standard hoc files in the model folder and copy them if not existing
if ~exist(fullfile(params.path,'lib_genroutines'),'file')
    mkdir(params.path,'lib_genroutines')
    display('non-existent folder lib_genroutines created')
end
if ~exist(fullfile(params.path,'lib_genroutines/fixnseg.hoc'),'file')
    copyfile(fullfile(t2npath,'fixnseg.hoc'),fullfile(params.path,'lib_genroutines/fixnseg.hoc'))
    display('fixnseg.hoc copied to model folder')
end
if ~exist(fullfile(params.path,'lib_genroutines/genroutines.hoc'),'file')
    copyfile(fullfile(t2npath,'genroutines.hoc'),fullfile(params.path,'lib_genroutines/genroutines.hoc'))
    display('genroutines.hoc copied to model folder')
end
if ~exist(fullfile(params.path,'lib_genroutines/pasroutines.hoc'),'file')
    copyfile(fullfile(t2npath,'pasroutines.hoc'),fullfile(params.path,'lib_genroutines/pasroutines.hoc'))
    display('pasroutines.hoc copied to model folder')
end


if params.cvode && isnumeric(params.dt)
    warning ('t2n:cvode', 'Dt is set but cvode is active. Dt will be ignored');
end

% create the local and server exchange folder
if exist(exchfolder,'dir') == 0
    mkdir(exchfolder);
end
if strfind(options,'-cl')
    [params.server.connect,outp] = sshfrommatlabissue(params.server.connect,sprintf('mkdir %s',params.server.modelfolder));  % check if mainfolder exists
    [params.server.connect,outp] = sshfrommatlabissue(params.server.connect,sprintf('rm -rf %s',nrn_exchfolder));
    [params.server.connect,outp] = sshfrommatlabissue(params.server.connect,sprintf('mkdir %s',nrn_exchfolder));
end


% check input
[ thesetrees,usestreesof ] = t2n_checkinput(neuron,tree);

%
for t = 1:numel(tree)
    if isfield(tree{t},'artificial') && ~isfield(tree{t},'NID')
        tree{t}.NID = strcat('cell_',tree{t}.artificial);           % artificial cells only need one morph hoc file which is named cell_ + the name of the artificial cell..
    end
end
NIDs = unique(cellfun(@(x) x.NID,tree,'UniformOutput',0));  % improves speed if many same cells are used
if ~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(morphfolder,strcat(x,'.hoc')),'file'),NIDs))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet or hoc file is missing! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        %         ind = cellfun(@(x) ~isfield(x,'NID'),tree); % find indices to not transformed trees
        ind = cellfun(@(x) ~isfield(x,'NID'),tree) ;
        if ~all(ind)
            ind(~ind) = ~cellfun(@(x) exist(fullfile(morphfolder,strcat(x.NID,'.hoc')),'file'),tree(~ind));
        end
        
        tree(ind) = t2n_writetrees(params,tree(ind),'',options);
    else
        error('T2N aborted');
    end
end
origminterf = cell(numel(tree),1);
for t = 1:numel(tree)
    if ~isfield(tree{t},'artificial')
        origminterf{t} = load(fullfile(morphfolder,sprintf('%s_minterf.dat',tree{t}.NID)));
    end
end


%% start writing hoc file
spines_flag = false(numel(neuron),1);
minterf = cell(numel(tree),1);
for n = 1:numel(neuron)
    x = getref(n,neuron,'mech');
    for t = 1:numel(neuron{x}.mech)   % check if any region of any tree has the spines mechanism
        if ~isempty(neuron{x}.mech{t})
            fields = fieldnames(neuron{x}.mech{t});
            for f = 1:numel(fields)
                if any(strcmpi(fieldnames(neuron{x}.mech{t}.(fields{f})),'spines'))
                    spines_flag(n) = true;
                    break
                end
            end
            if spines_flag(n)
                break
            end
        end
    end
    
    if ~isfield(neuron{n},'custom')
        neuron{n}.custom = {};
    end
    if strfind(options,'-d')
        tim = tic;
    end
    for tt = 1:numel(tree(thesetrees{n}))
        if ~isfield(tree{thesetrees{n}(tt)},'artificial')
            x = getref(n,neuron,'mech');
            minterf{thesetrees{n}(tt)} = make_nseg(tree{thesetrees{n}(tt)},origminterf{thesetrees{n}(tt)},params,neuron{x}.mech{thesetrees{n}(tt)});
        end
    end
    if ~isfield(params,'access')
        access = [find(~cellfun(@(y) isfield(y,'artificial'),tree(thesetrees{n})),1,'first'), 1];      % std accessing first non-artificial tree at node 1
    else
        access = params.access;
    end
    thisfolder = sprintf('sim%d',n);
    
    if exist(fullfile(exchfolder,thisfolder),'dir') == 0
        mkdir(fullfile(exchfolder,thisfolder));
    end
    if exist(fullfile(exchfolder,thisfolder,'iamrunning'),'file')
        answer = questdlg(sprintf('Error!\n%s seems to be run by another Matlab instance!\nOverwriting might cause errorneous output!\nIf you are sure that there is no simulation running, we can continue and overwrite. Are you sure? ',fullfile(exchfolder,thisfolder)),'Overwrite unfinished simulation','Yes to all','Yes','No (Cancel)','No (Cancel)');
        switch answer
            case 'Yes'
                % iamrunning file is kept and script goes on...
            case 'Yes to all'  % delete all iamrunning files in that exchfolder except from the current simulation
                folders = dir(exchfolder);
                for f = 1:numel(folders) % ignore first two as these are . and ..
                    if ~isempty(strfind(folders(f).name,'sim')) && ~strcmp(folders(f).name,thisfolder) && exist(fullfile(exchfolder,folders(f).name,'iamrunning'),'file')
                        delete(fullfile(exchfolder,folders(f).name,'iamrunning'))
                    end
                end
            otherwise
                error('T2N aborted')
        end
    else
        ofile = fopen(fullfile(exchfolder,thisfolder,'iamrunning') ,'wt');   %open morph hoc file in write modus
        fclose(ofile);
    end
    % delete the readyflag and log files if they exist
    if exist(fullfile(exchfolder,thisfolder,'readyflag'),'file')
        delete(fullfile(exchfolder,thisfolder,'readyflag'))
    end
    if exist(fullfile(exchfolder,thisfolder,'ErrorLogFile.txt'),'file')
        delete(fullfile(exchfolder,thisfolder,'ErrorLogFile.txt'))
    end
    if exist(fullfile(exchfolder,thisfolder,'NeuronLogFile.txt'),'file')
        delete(fullfile(exchfolder,thisfolder,'NeuronLogFile.txt'))
    end
    if ~isempty(strfind(options,'-cl'))
        [params.server.connect,outp] = sshfrommatlabissue(params.server.connect,sprintf('mkdir %s/%s',nrn_exchfolder,thisfolder));
    end
    
    %% write interface hoc
    
    nfile = fopen(fullfile(exchfolder,thisfolder,interf_file) ,'wt');   %open resulting hoc file in write modus
    
    fprintf(nfile,'// ***** This is a NEURON hoc file automatically created by the Matlab-NEURON interface. *****\n');
    fprintf(nfile,'// ***** Copyright by Marcel Beining and Johannes Kasper, Clinical Neuroanatomy, Goethe University Frankfurt*****\n\n');
    %initialize variables in NEURON
    fprintf(nfile,'// General variables: i, CELLINDEX, debug_mode, accuracy\n\n');
    
    fprintf(nfile,'// ***** Initialize Variables *****\n');
    fprintf(nfile,'strdef tmpstr // temporary string object\nobjref f\n');
    fprintf(nfile,'objref nil,cvode,strf,tvec,cell,cellList,pp,ppList,con,conList,nilcon,nilconList,rec,recList,rect,rectList,playt,playtList,play,playList,APCrec,APCrecList,APC,APCList,APCcon,APCconList,thissec,thisseg,thisval,maxRa,maxcm \n cellList = new List() // comprises all instances of cell templates, also artificial cells\n ppList = new List() // comprises all Point Processes of any cell\n conList = new List() // comprises all NetCon objects\n recList = new List() //comprises all recording vectors\n rectList = new List() //comprises all time vectors of recordings\n playtList = new List() //comprises all time vectors for play objects\n playList = new List() //comprises all vectors played into an object\n APCList = new List() //comprises all APC objects\n APCrecList = new List() //comprises all APC recording vectors\n nilconList = new List() //comprises all NULL object NetCons\n cvode = new CVode() //the Cvode object\n thissec = new Vector() //for reading range variables\n thisseg = new Vector() //for reading range variables\n thisval = new Vector() //for reading range variables\n\n');% maxRa = new Vector() //for reading range variables\n maxcm = new Vector() //for reading range variables\n\n');%[',numel(tree),']\n'  ;
    fprintf(nfile,sprintf('\nchdir("%s") // change directory to main simulation folder \n',nrn_path));
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Define some basic parameters *****\n');
    fprintf(nfile,sprintf('debug_mode = %d\n',debug) );
    if isfield(params,'accuracy')
        fprintf(nfile,sprintf('accuracy = %d\n',params.accuracy) );
    else
        fprintf(nfile,'accuracy = 0\n' );
    end
    fprintf(nfile,'strf = new StringFunctions()\n');
    if params.cvode
        fprintf(nfile,'cvode.active(1)\n');
        %         fprintf(nfile,'cvode.minstep(0.0000001)\n');    % less does not make any sense
        if params.use_local_dt
            fprintf(nfile,'io = cvode.use_local_dt(1)\n');
        end
    else
        fprintf(nfile,'io = cvode.active(0)\n');
        fprintf(nfile,sprintf('tvec = new Vector()\ntvec = tvec.indgen(%f,%f,%f)\n',params.tstart,params.tstop,params.dt));
        fprintf(nfile,'f = new File()\n');      %create a new filehandle
        fprintf(nfile,sprintf('io = f.wopen("%s//%s//tvec.dat")\n',params.exchfolder,thisfolder)  );  % open file for this time vector with write perm.
        fprintf(nfile,sprintf('io = tvec.printf(f,"%%%%-20.10g\\\\n")\n') );%"%%%%-20.10g")\n', c ) );    % print the data of the vector into the file
        fprintf(nfile,'io = f.close()\n');
    end
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Load standard libraries *****\n');
    if isfield(params,'nrnmech')
        if iscell(params.nrnmech)
            for c = 1:numel(params.nrnmech)
                if ~exist(fullfile(params.path,'lib_mech',params.nrnmech{c}),'file')
                    error('File %s is not existent in folder lib_mech',params.nrnmech{c})
                end
                fprintf(nfile,sprintf('io = nrn_load_dll("lib_mech/%s")\n',params.nrnmech{c}));
            end
        else
            if ~exist(fullfile(params.path,'lib_mech',params.nrnmech),'file')
                error('File %s is not existent in folder lib_mech',params.nrnmech)
            end
            fprintf(nfile,sprintf('io = nrn_load_dll("lib_mech/%s")\n',params.nrnmech));
        end
    else
        if ~exist(fullfile(params.path,'lib_mech','nrnmech.dll'),'file') || ~isempty(strfind(options,'-m'))  % check for existent file, otherwise compile dll
            nrn_installfolder = regexprep(fileparts(fileparts(params.neuronpath)),'\\','/');
            tstr = sprintf('cd "%s" && %s/mingw/bin/sh "%s/mknrndll.sh" %s',[nrn_path,'/lib_mech'],nrn_installfolder, regexprep(t2npath,'\\','/'), ['/',regexprep(nrn_installfolder,':','')]);
            [~,cmdout] = system(tstr);
            if isempty(strfind(cmdout,'nrnmech.dll was built successfully'))
                error('File nrnmech.dll was not found in lib_mech and compiling it with mknrndll failed! Check your mod files and run mknrndll manually')
            else
                display('nrnmech.dll compiled from mod files in folder lib_mech') 
            end
            rename_nrnmech()  % delete the o and c files
        end
        fprintf(nfile,'nrn_load_dll("lib_mech/nrnmech.dll")\n');
    end
    
    if params.openNeuron
        fprintf(nfile,'io = load_file("nrngui.hoc")\n');     % load the NEURON GUI
    else
        fprintf(nfile,'io = load_file("stdgui.hoc")\n');     % ony load other standard procedures
    end
    fprintf(nfile, sprintf('io = xopen("lib_genroutines/fixnseg.hoc")\n') );
    fprintf(nfile, sprintf('io = xopen("lib_genroutines/genroutines.hoc")\n') );
    fprintf(nfile, sprintf('io = xopen("lib_genroutines/pasroutines.hoc")\n') );
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Load custom libraries *****\n');
    if ~isempty(neuron{n}.custom)
        for c = 1:size(neuron{n}.custom,1)
            if strcmpi(neuron{n}.custom{c,2},'start')
                if strcmp(neuron{n}.custom{c,1}(end-4:end),'.hoc')   %check for hoc ending
                    if exist(fullfile(nrn_path,'lib_custom',neuron{n}.custom{c,1}),'file')
                        fprintf(nfile,sprintf('io = load_file("lib_custom/%s")\n',neuron{n}.custom{c,1}));
                    else
                        display(sprintf('File "%s" does not exist',neuron{n}.custom{c,1}))
                    end
                else
                    fprintf(nfile,neuron{n}.custom{c,1});  % add string as custom neuron code
                end
            end
        end
    end
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Load cell morphologies and create artificial cells *****\n');
    fprintf(nfile,sprintf('io = xopen("%s/%s/init_cells.hoc")\n',nrn_exchfolder,sprintf('sim%d',usestreesof(n)))); % das passt so!
    
    fprintf(nfile,sprintf('\n\nchdir("%s/%s") // change directory to folder of simulation #%d \n',params.exchfolder,thisfolder,n));
    fprintf(nfile,'\n\n');
    
    x = getref(n,neuron,'mech');
    if ~isnan(x)
        fprintf(nfile,'// ***** Load mechanisms and adjust nseg *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_mech.hoc")\n',nrn_exchfolder,sprintf('sim%d',x)) );
        else
            fprintf(nfile,'io = xopen("init_mech.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    
    x = getref(n,neuron,'pp');
    if ~isnan(x)
        fprintf(nfile,'// ***** Place Point Processes *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_pp.hoc")\n',nrn_exchfolder,sprintf('sim%d',x)) );
        else
            fprintf(nfile,'io = xopen("init_pp.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    x = getref(n,neuron,'con');
    if ~isnan(x)
        fprintf(nfile,'// ***** Define Connections *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_con.hoc")\n',nrn_exchfolder,sprintf('sim%d',x)) );
        else
            fprintf(nfile,'io = xopen("init_con.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    
    
    x = getref(n,neuron,'record');
    x2 = getref(n,neuron,'APCount');
    if ~isnan(x) || ~isnan(x2)
        fprintf(nfile,'// ***** Define recording sites *****\n');
        if x~=n && x2~=n && x==x2  % if both reference to the same other sim, use this sim
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_rec.hoc")\n',nrn_exchfolder,sprintf('sim%d',x)) );
        else  % else write an own file
            fprintf(nfile,'io = xopen("init_rec.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    
    x = getref(n,neuron,'play');
    if ~isnan(x)
        fprintf(nfile,'// ***** Define vector play sites *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_play.hoc")\n',nrn_exchfolder,sprintf('sim%d',x)) );
        else
            fprintf(nfile,'io = xopen("init_play.hoc")\n' );
        end
    end
    
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Last settings *****\n');
    if isfield(params,'celsius')
        if params.q10 == 1
            fprintf(nfile,'\n\nobjref q10\nq10 = new Temperature()\n' ) ;
            fprintf(nfile,sprintf('io = q10.correct(%g)\n\n',params.celsius) ) ;
        else
            fprintf(nfile,sprintf('celsius = %g\n\n',params.celsius) ) ;
            
        end
    end
    if spines_flag(n)
        fprintf(nfile,'addsurf_spines()\n');
    end
    fprintf(nfile,sprintf('tstart = %f\n',params.tstart));   %set params.tstart
    fprintf(nfile,sprintf('tstop = %f + %f //advances one more step due to roundoff errors for high tstops\n',params.tstop,params.dt));   %set params.tstop
    fprintf(nfile,sprintf('dt = %f\n',params.dt));         % set params.dt
    fprintf(nfile,sprintf('steps_per_ms = %f\n',1/params.dt));         % set steps per ms to avois changing dt on reinit
    if isfield(params,'v_init')
        fprintf(nfile,sprintf('v_init = %f\n',params.v_init));
    end
    fprintf(nfile,sprintf('prerun = %d\n',params.prerun));
    if numel(access) > 1 % if there is any non-artificial cell defined
        fprintf(nfile,sprintf('access cellList.o(%d).allregobj.o(%d).sec\n',access(1)-1,minterf{thesetrees{n}(access(1))}(access(2),2)) );
    end
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Include prerun or standard run replacing custom code *****\n');
    if ~isempty(neuron{n}.custom)
        for c = 1:size(neuron{n}.custom,1)
            if strcmpi(neuron{n}.custom{c,2},'mid')
                if strcmp(neuron{n}.custom{c,1}(end-4:end),'.hoc')   %check for hoc ending
                    if exist(fullfile(nrn_path,'lib_custom',neuron{n}.custom{c,1}),'file')
                        fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,neuron{n}.custom{c,1}));
                    else
                        display(sprintf('File "%s" does not exist',neuron{n}.custom{c,1}))
                    end
                else
                    fprintf(nfile,neuron{n}.custom{c,1});  % add string as custom neuron code
                end
            end
        end
    end
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Run NEURON *****\n');
    
    if ~params.skiprun
        fprintf(nfile,'init()\n');  % this needs to be modified later since v_init might be restarted
        fprintf(nfile,'run()\n');         % directly run the simulation
    else
        fprintf(nfile,'// Run is skipped due to custom code\n');
    end
    
    fprintf(nfile,'\n\n');
    fprintf(nfile,'io = xopen("save_rec.hoc")\n' );
    
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Include finishing custom code *****\n');
    if ~isempty(neuron{n}.custom)
        for c = 1:size(neuron{n}.custom,1)
            if strcmpi(neuron{n}.custom{c,2},'end')
                if strcmp(neuron{n}.custom{c,1}(end-4:end),'.hoc')   %check for hoc ending
                    if exist(fullfile(nrn_path,'lib_custom',neuron{n}.custom{c,1}),'file')
                        fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,neuron{n}.custom{c,1}));
                    else
                        display(sprintf('File "%s" does not exist',neuron{n}.custom{c,1}))
                    end
                else
                    fprintf(nfile,neuron{n}.custom{c,1});  % add string as custom neuron code
                end
            end
        end
    end
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Make Matlab notice end of simulation *****\n');
    fprintf(nfile,'f = new File()\n');       %create a new filehandle
    fprintf(nfile,'io = f.wopen("readyflag")\n' );       % create the readyflag file
    fprintf(nfile,'io = f.close()\n');   % close the filehandle
    if ~params.openNeuron
        fprintf(nfile,'quit()\n');  % exit NEURON if it was defined so in the parameters
    end
    
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// *-*-*-*-* END *-*-*-*-*\n');
    
    fclose(nfile);
    
    
    %% write init_cells.hoc
    
    if usestreesof(n) == n  % write only if morphologies are not referenced to other sim init_cell.hoc
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_cells.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Load cell morphology templates and create artificial cells *****\n');
        templates = cell(0);
        for tt = 1:numel(thesetrees{n})
            % load templates generated by neuron_template_tree, create one
            % instance of them and add them to the cellList
            if ~any(strcmp(templates,tree{thesetrees{n}(tt)}.NID))
                fprintf(ofile,sprintf('io = xopen("%s//%s.hoc")\n',params.morphfolder,tree{thesetrees{n}(tt)}.NID) );
                templates = cat(1,templates,tree{thesetrees{n}(tt)}.NID);
            end
            fprintf(ofile, sprintf('cell = new %s()\n', tree{thesetrees{n}(tt)}.NID) );
            if isfield( tree{thesetrees{n}(tt)},'params')
                fields = fieldnames( tree{thesetrees{n}(tt)}.params);
                for f = 1:numel(fields)
                    if ischar(tree{thesetrees{n}(tt)}.params.(fields{f})) && regexpi(tree{thesetrees{n}(tt)}.params.(fields{f}),'^(.*)$')  % check if this is a class/value pair, then use the () instead of =
                        fprintf(ofile, sprintf('cell.cell.%s%s\n',fields{f}, tree{thesetrees{n}(tt)}.params.(fields{f})));
                    else
                        fprintf(ofile, sprintf('cell.cell.%s = %g\n',fields{f}, tree{thesetrees{n}(tt)}.params.(fields{f})));
                    end
                end
            end
            
            fprintf(ofile, 'io = cellList.append(cell)\n');
            
        end
        fprintf(ofile, 'objref cell\n');
        
        fprintf(ofile,'\n\n');
        
        fclose(ofile);
    elseif exist(fullfile(exchfolder,thisfolder,'init_cells.hoc'),'file')
        delete(fullfile(exchfolder,thisfolder,'init_cells.hoc'));
    end
    
    %% write init_mech.hoc
    
    if getref(n,neuron,'mech') == n     %rewrite only if mechanism is not taken from previous sim
        rangestr = '';
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_mech.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Insert mechanisms *****\n');
        if isfield(neuron{n},'mech')
            flag_strion_all = false(numel(thesetrees{n}),1);
            strion_all = cell(numel(thesetrees{n}),1);
            flag_strnseg_all = false(numel(thesetrees{n}),1);
            strnseg_all = cell(numel(thesetrees{n}),1);
            strion_reg = cell(numel(thesetrees{n}),1);
            flag_strion_reg = cell(numel(thesetrees{n}),1);
            strnseg_reg = cell(numel(thesetrees{n}),1);
            flag_strnseg_reg = cell(numel(thesetrees{n}),1);
            for tt = 1:numel(thesetrees{n})
                t = thesetrees{n}(tt);
                if numel(neuron{n}.mech) >= t && ~isempty(neuron{n}.mech{t})   && ~isfield(tree{thesetrees{n}(tt)},'artificial')    % if a mechanism is defined for this tree
                    if isstruct(neuron{n}.mech{t})          % input must be a structure
                        fields = fieldnames(neuron{n}.mech{t});
                    else
                        continue
                    end
                    
                    if any(strcmpi(fields,'all'))
                        str = sprintf('forsec cellList.o(%d).allreg {\n',tt-1);   %neuron:go through this region
                        strion_all{tt} = sprintf('forsec cellList.o(%d).allreg {\n',tt-1);   %neuron:go through this region
                        
                        strnseg_all{tt} = sprintf('forsec cellList.o(%d).allreg {\n',tt-1);   %neuron:go through this region
                        
                        mechs = fieldnames(neuron{n}.mech{t}.all);                % mechanism names are the fieldnames in the structure
                        if any(strcmp(mechs,'nseg'))
                            mechs = setdiff(mechs,'nseg');
                            strnseg_all{tt} = sprintf('%snseg = %d\n',strnseg_all{tt},neuron{n}.mech{t}.all.nseg);   %neuron: define values
                            flag_strnseg_all(tt) = true;
                        end
                        for m = 1:numel(mechs)      % loop through mechanisms
                            str = sprintf('%sinsert %s\n',str,mechs{m});        % neuron:insert this mechanism
                            if ~isempty(neuron{n}.mech{t}.all.(mechs{m}))
                                mechpar = fieldnames(neuron{n}.mech{t}.all.(mechs{m}));
                                for p = 1:numel(mechpar)  % loop through mechanism parameters
                                    if strcmpi(mechpar{p},'cm') || strcmpi(mechpar{p},'Ra') || (~isempty(strfind(mechs{m},'_ion')) &&  (numel(mechpar{p}) <= strfind(mechs{m},'_ion') || (numel(mechpar{p}) > strfind(mechs{m},'_ion') && ~strcmp(mechpar{p}(strfind(mechs{m},'_ion')+1),'0'))))       %if mechanism is an ion or passive cm/Ra, leave out mechansim suffix
                                        if ~isempty(strfind(mechs{m},'_ion')) && strcmpi(mechpar{p},'style')
                                            if numel(neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p})) ~= 5
                                                for nn = 1:numel(neuron)
                                                    delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                end
                                                error(sprintf('Error! Style specification of ion "%s" should be 5 numbers (see NEURON or t2n documentation)',mechs{m}))

%                                                 return
                                            end
                                            strion_all{tt} = sprintf('%sion_style("%s",%d,%d,%d,%d,%d)\n',strion_all{tt},mechs{m},neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p}));   %neuron: define values
                                            flag_strion_all(tt) = true;
                                        else
                                            str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p}));   %neuron: define values
                                        end
                                    else
                                        str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p}));   %neuron: define values
                                    end
                                end
                            end
                        end
                        fprintf(ofile,sprintf('%s}\n\n',str));
                    end
                    
                    if isfield(tree{thesetrees{n}(tt)},'R')
                        uR = unique(tree{thesetrees{n}(tt)}.R); % Region indices that exist in tree
                        if ~isempty(intersect(tree{thesetrees{n}(tt)}.rnames(uR),fields)) %isstruct(neuron{n}.mech{t}.(fields{1}))  %check if mechanism was defined dependent on existent region
                            regs = fields;  %if yes (some of) the input are the regions
                            regs = intersect(tree{thesetrees{n}(tt)}.rnames(uR),regs);  % only use those region names which are existent in tree
                            strion_reg{tt} = cell(numel(regs),1);
                            flag_strion_reg{tt} = false(numel(regs),1);
                            strnseg_reg{tt} = cell(numel(regs),1);
                            flag_strnseg_reg{tt} = false(numel(regs),1);
                            for r = 1 : numel(regs)
                                str = sprintf('forsec cellList.o(%d).reg%s {\n',tt-1,regs{r});   %neuron:go through this region
                                strnseg_reg{tt}{r} = sprintf('forsec cellList.o(%d).reg%s {\n',tt-1,regs{r});   %neuron:go through this region
                                mechs = fieldnames(neuron{n}.mech{t}.(regs{r}));                % mechanism names are the fieldnames in the structure
                                if any(strcmp(mechs,'nseg'))
                                    mechs = setdiff(mechs,'nseg');
                                    strnseg_reg{tt}{r} = sprintf('%snseg = %d\n',strnseg_reg{tt}{r},neuron{n}.mech{t}.(regs{r}).nseg);   %neuron: define values
                                    flag_strnseg_reg{tt}(r) = true;
                                end
                                for m = 1:numel(mechs)      % loop through mechanisms
                                    str = sprintf('%sinsert %s\n',str,mechs{m});        % neuron:insert this mechanism
                                    
                                    if ~isempty(neuron{n}.mech{t}.(regs{r}).(mechs{m}))
                                        mechpar = fieldnames(neuron{n}.mech{t}.(regs{r}).(mechs{m}));
                                        for p = 1:numel(mechpar)  % loop through mechanism parameters
                                            if strcmpi(mechpar{p},'cm') || strcmpi(mechpar{p},'Ra') || (~isempty(strfind(mechs{m},'_ion')) &&  (numel(mechpar{p}) <= strfind(mechs{m},'_ion') || (numel(mechpar{p}) > strfind(mechs{m},'_ion') && ~strcmp(mechpar{p}(strfind(mechs{m},'_ion')+1),'0'))))       %if mechanism is an ion or passive cm/Ra, leave out mechansim suffix
                                                if ~isempty(strfind(mechs{m},'_ion')) && strcmpi(mechpar{p},'style')
                                                    if numel(neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p})) ~= 5
                                                        for nn = 1:numel(neuron)
                                                            delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                        end
                                                        error(sprintf('Error! Style specification of ion "%s" should be 5 numbers (see NEURON or t2n documentation)',mechs{m}))

%                                                         return
                                                    end
                                                    strion_reg{tt}{r} = sprintf('%sion_style("%s",%d,%d,%d,%d,%d)\n',strion_reg{tt}{r},mechs{m},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                                    flag_strion_reg{tt}(r) = true;
                                                else
                                                    str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                                end
                                                
                                            else
                                                if numel(neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p})) == 1
                                                    str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                                else
                                                    for nn = 1:numel(neuron)
                                                        delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                    end
                                                    error(sprintf('Parameter %s of mechanism %s in region %s has more than one value, please check.',mechpar{p},mechs{m},regs{r}))
%                                                     return
                                                end
                                            end
                                        end
                                    end
                                end
                                fprintf(ofile,sprintf('%s}\n\n',str));
                            end
                        end
                    end
                    if any(strcmpi(fields,'range'))
                        if ~isfield(tree{thesetrees{n}(tt)},'artificial')
                            %                             str = '';
                            [~, ia] = unique(minterf{thesetrees{n}(tt)}(:,[2,4]),'rows','stable');  % find real segments in neuron simulation
                            ia = ia(~isnan(minterf{thesetrees{n}(tt)}(ia,4))); % remove start nodes of a segment (cause their value belongs to segment -1)
                            ia(numel(ia)+1) = size(minterf{thesetrees{n}(tt)},1)+1;   % add one entry
                            
                            mechs = fieldnames(neuron{n}.mech{t}.range);
                            for m = 1:numel(mechs)
                                vars = fieldnames(neuron{n}.mech{t}.range.(mechs{m}));
                                %                                 allvals = zeros(3,0);
                                %                                 thesevars = '';
                                for r = 1:numel(vars)
                                    if numel(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})) == numel(tree{thesetrees{n}(tt)}.X)
                                        allvals = zeros(3,0);
                                        for in = 1:numel(ia)-1
                                            thisval = nanmean(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})(minterf{thesetrees{n}(tt)}(ia(in),1):minterf{thesetrees{n}(tt)}(ia(in+1)-1,1))); % the value is the mean of all tree nodes which are simulated by this segment, if first node is start of section, ignore this one, since it belongs to old region
                                            
                                            if ~isnan(thisval)
                                                allvals = cat(2,allvals,[minterf{thesetrees{n}(tt)}(ia(in),[2,4]),thisval]');
                                            end
                                        end
                                        
                                        %                                         thesevars = sprintf('%s"%s_%s",',thesevars,vars{r},mechs{m});
                                        secname = sprintf('range_%s_%s_%s_sec.dat',tree{thesetrees{n}(tt)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,secname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(1,:));
                                        fclose(f);
                                        segname = sprintf('range_%s_%s_%s_seg.dat',tree{thesetrees{n}(tt)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,segname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(2,:));
                                        fclose(f);
                                        valname = sprintf('range_%s_%s_%s_val.dat',tree{thesetrees{n}(tt)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,valname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(3,:));
                                        fclose(f);
                                        if any(strcmp({'cm','Ra'},vars{r}))  % if variable is cm or Ra, do not write _"mech"  behind it
                                            rangestr = sprintf('%sset_range(%d,"%s","%s","%s","%s")\n',rangestr,tt-1,secname,segname,valname,vars{r});
                                        else
                                            rangestr = sprintf('%sset_range(%d,"%s","%s","%s","%s_%s")\n',rangestr,tt-1,secname,segname,valname,vars{r},mechs{m});
                                        end
                                    else
                                        for nn = 1:numel(neuron)
                                            delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                        end
                                        error('Range variable definition should be a vector with same number of elements as tree has nodes')
%                                         return
                                    end
                                end
                            end
                            
                        else
                            for nn = 1:numel(neuron)
                                delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                            end
                            error('Setting range variables for artificial cells is invalid')
                        end
                    end
                end
            end
            fprintf(ofile,'\n\n');
            
            if any(cellfun(@any,flag_strion_reg)) || any(flag_strion_all)
                fprintf(ofile,'// ***** Now add specific ion styles *****\n');
                if any(flag_strion_all)
                    for tt = 1:numel(thesetrees{n})
                        if flag_strion_all(tt)
                            fprintf(ofile,sprintf('%s}\n\n',strion_all{tt}));
                        end
                    end
                end
                if any(cellfun(@any,flag_strion_reg))
                    for tt = 1:numel(thesetrees{n})
                        for r = 1:numel(flag_strion_reg{tt})
                            if flag_strion_reg{tt}(r)
                                fprintf(ofile,sprintf('%s}\n\n',strion_reg{tt}{r}));
                            end
                        end
                    end
                end
            end
            
            fprintf(ofile,'// ***** Define nseg for all cells *****\n');
            fprintf(ofile, 'proc make_nseg() {\n');
            fprintf(ofile, 'for CELLINDEX = 0, cellList.count -1 {\n');
            fprintf(ofile, 'if (cellList.o(CELLINDEX).is_artificial == 0) {\n');
            if isfield(params,'nseg') && isnumeric(params.nseg)
                fprintf(ofile, 'forsec cellList.o(CELLINDEX).allreg {\n');
                fprintf(ofile, sprintf('nseg = %f\n}\n}\n}\n\n',round(params.nseg)) );
                if rem(round(params.nseg),2) == 0
                    warning('nseg is not odd! Please reconsider nseg');
                end
            elseif isfield(params,'nseg') && strcmpi(params.nseg,'dlambda')
                fprintf(ofile, 'geom_nseg()\n}\n}\n\n');
            elseif isfield(params,'nseg') && ~isempty(strfind(params.nseg,'ach'))
                each = cell2mat(textscan(params.nseg,'%*s %d')); % get number
                fprintf(ofile, 'forsec cellList.o(CELLINDEX).allreg {\n');
                fprintf(ofile, sprintf('n = L/%d\nnseg = n+1\n}\n}\n}\n\n',each) );
            else
                fprintf(ofile, '// No nseg specified!!!\n}\n}\n\n');
                warning('nseg has not been specified in params.nseg (correctly?)!')
            end
            
            if any(cellfun(@any,flag_strnseg_reg)) || any(flag_strnseg_all)
                fprintf(ofile,'// ***** Add specific nseg definitions *****\n');
                if any(flag_strnseg_all)
                    for tt = 1:numel(thesetrees{n})
                        fprintf(ofile,sprintf('%s}\n\n',strnseg_all{tt}));
                    end
                end
                if any(cellfun(@any,flag_strnseg_reg))
                    for tt = 1:numel(thesetrees{n})
                        for r = 1:numel(flag_strnseg_reg{tt})
                            if flag_strnseg_reg{tt}(r)
                                fprintf(ofile,sprintf('%s}\n\n',strnseg_reg{tt}{r}));
                            end
                        end
                    end
                end
            end
            fprintf(ofile,'}\n\n');
            
            
            fprintf(ofile,'// ***** Now adjust number of segments *****\n');
            fprintf(ofile,'make_nseg()\n');
            if ~isempty(rangestr)
                fprintf(ofile,'\n\n');
                fprintf(ofile,'// ***** Set specified range variables *****\n');
                fprintf(ofile,rangestr);
                fprintf(ofile,'\n\n');
            end
        end
        fclose(ofile);          %close file
    elseif exist(fullfile(exchfolder,thisfolder,'init_mech.hoc'),'file')
        delete(fullfile(exchfolder,thisfolder,'init_mech.hoc'));
    end
    
    %% write init_pp.hoc
    
    if getref(n,neuron,'pp') == n     %rewrite only if PP def is not taken from previous sim
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_pp.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Place synapses, electrodes or other point processes *****\n');
        if isfield(neuron{n},'pp')
            count = 0;%zeros(numel(thesetrees{n}),1);
            for tt = 1:numel(thesetrees{n})
                t = thesetrees{n}(tt);
                if numel(neuron{n}.pp) >= t && ~isempty(neuron{n}.pp{t})   && ~isfield(tree{thesetrees{n}(tt)},'artificial')    % if point processes are defined for this tree
                    ppfield = fieldnames(neuron{n}.pp{t});
                    for f1 = 1:numel(ppfield)
                        %%%%
                        for n1 = 1:numel(neuron{n}.pp{t}.(ppfield{f1}))
                            node = neuron{n}.pp{t}.(ppfield{f1})(n1).node;
                            %                         if numel(node) < 5
                            %
                            %                         else
                            %                             secname = sprintf('pp_%s_%s_%s_sec.dat',ppfield{f1},tree{thesetrees{n}(tt)}.NID);
                            %                                         f = fopen(fullfile(exchfolder,thisfolder,secname) ,'Wt');
                            %                                         fprintf(f,'%g\n',allvals(1,:));
                            %                                         fclose(f);
                            %                                         segname = sprintf('range_%s_%s_%s_seg.dat',tree{thesetrees{n}(tt)}.NID,vars{r},mechs{m});
                            %                                         f = fopen(fullfile(exchfolder,thisfolder,segname) ,'Wt');
                            %                                         fprintf(f,'%g\n',allvals(2,:));
                            %                                         fclose(f);
                            %                         end
                            for in = 1:numel(node)
                                inode = find(minterf{thesetrees{n}(tt)}(:,1) == node(in),1,'first');    %find the index of the node in minterf
                                fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',tt-1,minterf{thesetrees{n}(tt)}(inode,2) ) );    % corresponding section of node
                                fprintf(ofile,sprintf('{pp = new %s(%f)\n',ppfield{f1},minterf{thesetrees{n}(tt)}(inode,3) ) );  % new pp
                                fields = setdiff(fieldnames(neuron{n}.pp{t}.(ppfield{f1})(n1)),{'node','id'});
                                
                                if any(strcmp(ppfield{f1},{'IClamp','SEClamp','SEClamp2','VClamp'})) && (any(strcmp(fields,'times')) || (any(strcmp(fields,'dur')) && (numel(neuron{n}.pp{t}.(ppfield{f1})(n1).dur) > 3 || neuron{n}.pp{t}.(ppfield{f1})(n1).del == 0)))  % check if field "times" exists or multiple durations are given or del is zero (last point can introduce a bug when cvode is active)
                                    if any(strcmp(fields,'times'))
                                        times = sprintf('%f,',neuron{n}.pp{t}.(ppfield{f1})(n1).times);
                                    else   %bugfix since seclamp can only use up to 3 duration specifications
                                        times = sprintf('%f,',[0 cumsum(neuron{n}.pp{t}.(ppfield{f1})(n1).dur(1:end-1))]);
                                    end
                                    amps = sprintf('%f,',neuron{n}.pp{t}.(ppfield{f1})(n1).amp);
                                    times = times(1:end-1); amps = amps(1:end-1); % delete last commas
                                    
                                    fprintf(ofile,'playt = new Vector()\n');
                                    fprintf(ofile,sprintf('playt = playt.append(%s)\n',times));
                                    fprintf(ofile,'play = new Vector()\n');
                                    fprintf(ofile,sprintf('play = play.append(%s)\n',amps));
                                    
                                    switch ppfield{f1}
                                        case 'IClamp'    % if SEClamp and VClamp dur and amp would be handled equally this could be simplified much more =/
                                            fprintf(ofile,'play.play(&pp.amp,playt)\n');
                                            fprintf(ofile,'pp.dur = 1e15\n');
                                            fprintf(ofile,'pp.del = -1e4\n');
                                        case 'VClamp'
                                            fprintf(ofile,'play.play(&pp.amp[0],playt)\n');
                                            fprintf(ofile,'pp.dur[0] = 1e15\n');
                                        case {'SEClamp','SEClamp2'}
                                            fprintf(ofile,'play.play(&pp.amp1,playt)\n');
                                            fprintf(ofile,'pp.dur1 = 1e15\n');
                                    end
                                    fields = setdiff(fields,{'times','amp','dur','del'});
                                    fprintf(ofile,'io = playtList.append(playt)\n');
                                    fprintf(ofile,'io = playList.append(play)\n');
                                    fprintf(ofile, 'objref play\n');
                                    fprintf(ofile, 'objref playt\n');
                                end
                                
                                for f2 =1:numel(fields)  % go through all parameter fields and declare them
                                    if any(strcmpi(fields{f2},{'dur','amp'})) && any(strcmp(ppfield{f1},{'SEClamp2','SEClamp','VClamp'}))   % for dur and amp, there are multiple values
                                        for ff = 1:numel(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}))
                                            switch ppfield{f1}
                                                case 'VClamp'
                                                    fprintf(ofile,sprintf('pp.%s[%d] = %f \n',fields{f2},ff-1,neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})(ff)));
                                                case {'SEClamp','SEClamp2'}
                                                    fprintf(ofile,sprintf('pp.%s%d = %f \n',fields{f2},ff,neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})(ff)) );
                                            end
                                        end
                                        
                                    else
                                        if numel(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})) > 1
                                            if numel(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})) == numel(node)
                                                if iscell(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})(in)) && ischar(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}){in}) && regexpi(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}){in},'^(.*)$')  % check if this is a class/value pair, then use the () instead of =
                                                    fprintf(ofile,sprintf('pp.%s%s \n', fields{f2},neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}){in} ));
                                                else
                                                    fprintf(ofile,sprintf('pp.%s = %f \n', fields{f2},neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})(in)) );
                                                end
                                            else
                                                for nn = 1:numel(neuron)
                                                    delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                end
                                                error(sprintf('Caution: "%s" vector of PP "%s" has wrong size!\n It has to be equal 1 or equal the number of nodes where the PP is inserted,',fields{f2},ppfield{f1}))
                                            end
                                        else
                                            if ischar(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})) && regexpi(neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}),'^(.*)$')  % check if this is a class/value pair, then use the () instead of =
                                                fprintf(ofile,sprintf('pp.%s%s \n', fields{f2},neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2}) ));
                                            else
                                                fprintf(ofile,sprintf('pp.%s = %f \n', fields{f2},neuron{n}.pp{t}.(ppfield{f1})(n1).(fields{f2})) );
                                            end
                                        end
                                    end
                                end
                                
                                fprintf(ofile,'}\n');
                                fprintf(ofile,'io = ppList.append(pp)\n' );  %append pp to ppList
                                neuron{n}.pp{t}.(ppfield{f1})(n1).id(in) = count;   % save id to pplist in Neuron (for find the correct object for recording later)
                                count = count +1; %ppnum(t) = ppnum(t) +1;  % counter up
                            end
                        end
                    end
                end
            end
            fprintf(ofile, 'objref pp\n');
        end
        fclose(ofile);
    end
    
    
    %% write init_con.hoc
    
    if getref(n,neuron,'con') == n     %rewrite only if connections are not taken from previous sim
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_con.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define Connections *****\n');
        if isfield(neuron{n},'con')
            % should have look like: {source(node or point process), what to
            % watch, target, threshold, delay, weight}
            for c = 1:numel(neuron{n}.con)
                str = cell(0);
                nodeflag = false;
                sourcefields = setdiff(fieldnames(neuron{n}.con(c).source),{'cell','watch'});
                
                cell_source = neuron{n}.con(c).source.cell;
                if isempty(sourcefields)   % probably an artificial cell...in that case "cell_source" can be a multi array, create a NetCon for each of these sources
                    for t = 1:numel(cell_source)
                        if ~isempty(cell_source(t)) && isfield(tree{cell_source(t)},'artificial')
                            str{t} = sprintf('con = new NetCon(cellList.o(%d).cell,',find(thesetrees{n}==cell_source(t))-1);
                            
                        else
                            str{t} = sprintf('con = new NetCon(nil,');
                        end
                    end
                else
                    if numel(cell_source) > 1
                        for nn = 1:numel(neuron)
                            delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                        end
                        error(sprintf('Error in connection %d of neuron instance %d! You must define a single cell ID if the source of a NetCon is a PointProcess or a section!',c,n))
                    end
                    if any(strcmp(sourcefields,'pp'))  % point process is the source
                        x = getref(n,neuron,'pp');
                        pp = neuron{n}.con(c).source.pp;
                        %                         [~,iid] = intersect(neuron{x}.pp{cell_source}.(pp).node,neuron{n}.con(c).source.node); % get reference to the node location of the PPs that should be connected
                        % that is not working if several pp are defined at
                        % that node. here comes the workaround
                        if isfield(neuron{n}.con(c).source,'ppg')  % check for an index to a PP subgroup
                            ppg = neuron{n}.con(c).source.ppg;
                        else
                            ppg = 1:numel(neuron{x}.pp{cell_source}.(pp));  % else take all PP subgroups
                        end
                        
                        upp = unique(cat(1,neuron{x}.pp{cell_source}.(pp)(ppg).node));  % unique pp nodes of the cell
                        if numel(upp) == 1 % otherwise hist would make as many bins as upp
                            cpp = numel(cat(1,neuron{x}.pp{cell_source}.(pp)(ppg).node));%1;
                        else
                            cpp =hist(cat(1,neuron{x}.pp{cell_source}.(pp)(ppg).node),upp); % number of pps at the same nodes
                        end
                        ucon = unique(neuron{n}.con(c).source.node); % unique connection nodes declared to the cell
                        if numel(ucon) == 1  % otherwise hist would make as many bins as ucon
                            ccon = numel(neuron{n}.con(c).source.node);
                        else
                            ccon =hist(neuron{n}.con(c).source.node,ucon); % number of connections to the same node
                        end
                        iid = cell(numel(neuron{n}.pp{cell_source}.(pp)(ppg)),1);  % initialize ids to pps
                        for uc = 1:numel(ucon)  % go trough all nodes that should be connected from
                            if any(ucon(uc) == upp)  % check if the pp exists there
                                for n1 = 1:numel(iid)
                                    ind = find(neuron{x}.pp{cell_source}.(pp)(ppg(n1)).node == ucon(uc));  % find all PPs at that node
                                    if ~isempty(ind)
                                        if cpp(ucon(uc) == upp) == ccon(uc)    % same number of PPs and connections, put them together, should be ok without warning
                                            iid{n1} = cat (1,iid{n1},ind);
                                        else
                                            for nn = 1:numel(neuron)
                                                delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                            end
                                            error(sprintf('Error cell %d. %d connections are declared to start from from %d %ss at node %d. Making a connection from PP at a node where multiple of these PPs exist is not allowed. Probably you messed something up',cell_source,ccon(uc),cpp(ucon(uc) == upp),pp,ucon(uc))) % give a warning if more connections were declared than PPs exist at that node
                                        end
                                    end
                                end
                                
                            else
                                display(sprintf('Warning cell %d. PP %s for connection does not exist at node %d',neuron{n}.con(c).target.cell,pp,ucon(uc)))
                            end
                        end
                        
                        for ii = 1:numel(iid)
                            str{ii} = sprintf('con = new NetCon(ppList.o(%d),',neuron{x}.pp{cell_source}.(pp).id(iid(ii)));
                        end
                    else   % a normal section is the source
                        node = neuron{n}.con(c).source.node;
                        for in = 1:numel(node)
                            inode = find(minterf{cell_source}(:,1) == node(in),1,'first');    %find the index of the node in minterf
                            if isfield(neuron{n}.con(c).source,'watch') && ischar(neuron{n}.con(c).source.watch)
                                str{in} = sprintf('cellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&%s(%f),',find(thesetrees{n}==cell_source)-1,minterf{thesetrees{n}(cell_source)}(inode,2),neuron{n}.con(c).source.watch,minterf{thesetrees{n}==cell_source}(inode,3));
                            else
                                str{in} = sprintf('cellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&v(%f),',find(thesetrees{n}==cell_source)-1,minterf{thesetrees{n}(cell_source)}(inode,2),minterf{thesetrees{n}==cell_source}(inode,3));
                            end
                        end
                        nodeflag = true;
                    end
                    
                end
                
                %%%
                
                targetfields = setdiff(fieldnames(neuron{n}.con(c).target),'cell');
                newstr = cell(0);
                count = 1;
                for it = 1:numel(neuron{n}.con(c).target)
                    cell_target = neuron{n}.con(c).target(it).cell;
                    if isempty(targetfields)   % probably an artificial cell...
                        for t1 = 1:numel(cell_source)
                            for t2 = 1:numel(cell_target)
                                if ~isempty(cell_target(t2)) && isfield(tree{cell_target(t2)},'artificial')
                                    newstr{count} = sprintf('%scellList.o(%d).cell',str{t1},find(thesetrees{n}==cell_target(t2))-1);
                                else
                                    newstr{count} = sprintf('%snil',str{t1});
                                end
                                count = count +1;
                            end
                        end
                    elseif any(strcmp(targetfields,'pp'))  % point process is the target
                        x = getref(n,neuron,'pp');
                        
                        pp = neuron{n}.con(c).target(it).pp;
                        %                         %                         neuron{x}.pp{t_target}.(pp).node
                        %                         [~,iid] = intersect(neuron{x}.pp{cell_target}.(pp).node,neuron{n}.con(c).target(it).node); % get reference to the node location of the PPs that should be connected
                        %                         intersect unfortunately fails if more than one PP
                        %                         is declared at the same node and con wants to
                        %                         target both. here comes the workaround
                        if isfield(neuron{n}.con(c).target(it),'ppg')
                            ppg = neuron{n}.con(c).target(it).ppg;
                        else
                            ppg = 1:numel(neuron{x}.pp{cell_target}.(pp));
                        end
                        upp = unique(cat(1,neuron{x}.pp{cell_target}.(pp)(ppg).node));  % unique pp nodes of the cell
                        if numel(upp) == 1 % otherwise hist would make as many bins as upp
                            cpp = numel(cat(1,neuron{x}.pp{cell_target}.(pp)(ppg).node));%1;
                        else
                            cpp =hist(cat(1,neuron{x}.pp{cell_target}.(pp)(ppg).node),upp); % number of pps at the same nodes
                        end
                        ucon = unique(neuron{n}.con(c).target(it).node); % unique connection nodes declared to the cell
                        if numel(ucon) == 1  % otherwise hist would make as many bins as ucon
                            ccon = numel(neuron{n}.con(c).target(it).node);
                        else
                            ccon =hist(neuron{n}.con(c).target(it).node,ucon); % number of connections to the same node
                        end
                        iid = cell(numel(neuron{n}.pp{cell_target}.(pp)(ppg)),1);  % initialize ids to pps
                        for uc = 1:numel(ucon)  % go trough all nodes that should be connected to
                            if any(ucon(uc) == upp)  % check if the pp exists there
%                                 
                                for n1 = 1:numel(iid)
                                    ind = find(neuron{x}.pp{cell_target}.(pp)(ppg(n1)).node == ucon(uc));  % find all PPs at that node
                                    if ~isempty(ind)
                                        if cpp(ucon(uc) == upp) < ccon(uc)  % less PPs exist than connections to node
                                            if numel(ind) == 1
                                                iid{n1} = cat (1,iid{n1},repmat(ind,ccon(uc),1));  % add as many PPs from that node to the id list as connections were declared (or as pps exist there)
                                                display(sprintf('Warning cell %d. More connections to same %s declared than %ss at that node. All connections target now that %s.',cell_target,pp,pp,pp)) % give a warning if more connections were declared than PPs exist at that node
                                            else
                                                for nn = 1:numel(neuron)
                                                    delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                end
                                                error(sprintf('Error cell %d. %d connections are declared to %d %ss at node %d. Probably you messed something up',cell_target,ccon(uc),cpp(ucon(uc) == upp),pp,ucon(uc))) % give a warning if more connections were declared than PPs exist at that node
                                            end
                                        elseif cpp(ucon(uc) == upp) > ccon(uc)  % more PPs exist than connections to node
                                            if isfield(neuron{n}.con(c).target(it),'ipp')
                                                iid{n1} = cat (1,iid{n1},ind(neuron{n}.con(c).target(it).ipp));
                                            else
                                                iid{n1} = cat (1,iid{n1},ind(1:min(cpp(ucon(uc) == upp),ccon(uc))));  % add as many PPs from that node to the id list as connections were declared (or as pps exist there)
                                                display(sprintf('Warning cell %d, node %d. Less connections to same %s declared than %ss at that node. Connections target now only the first %d %ss.',cell_target,ucon(uc),pp,pp,min(cpp(ucon(uc) == upp),ccon(uc)),pp)) % give a warning if more connections were declared than PPs exist at that node
                                            end
                                        else   % same number of PPs and connections, put them together, should be ok without warning
                                            iid{n1} = cat (1,iid{n1},ind);
                                        end
                                        
                                    end
                                end
                                
                            else
                                display(sprintf('Warning cell %d. PP %s for connection does not exist at node %d',cell_target,pp,ucon(uc)))
                            end
                        end
                        if numel(unique(cat(1,iid{:}))) ~= numel(cat(1,iid{:}))
                           display(sprintf('Warning cell %d. Connection #%d targets the PP %s at one or more nodes where several %s groups are defined! Connection is established to all of them. Use "neuron.con(x).target(y).ppg = z" to connect only to the zth group of PP %s.',neuron{n}.con(c).target.cell,c,pp,pp,pp))
                        end
                        for t1 = 1:numel(cell_source)
                            for n1 = 1:numel(iid)
                                for ii = 1:numel(iid{n1})
                                    newstr{count} = sprintf('%sppList.o(%d)',str{t1},neuron{x}.pp{cell_target}.(pp)(ppg(n1)).id(iid{n1}(ii)));
                                    count = count +1;
                                end
                            end
                        end
                        
                    else   % nothing...
                        warning('No target specified as connection')
                        for t1 = 1:numel(cell_source)
                            newstr{count} = sprintf('%snil',str{t1});
                            count = count + 1;
                        end
                    end
                end
                for s = 1:numel(newstr)
                    if isfield(neuron{n}.con(c),'threshold')
                        newstr{s} = sprintf('%s,%g', newstr{s},neuron{n}.con(c).threshold);
                    else
                        newstr{s} = sprintf('%s,10', newstr{s});
                    end
                    if isfield(neuron{n}.con(c),'delay')
                        newstr{s} = sprintf('%s,%g', newstr{s},neuron{n}.con(c).delay);
                    else
                        newstr{s} = sprintf('%s,1', newstr{s});
                    end
                    if isfield(neuron{n}.con(c),'weight')
                        newstr{s} = sprintf('%s,%g)\n', newstr{s},neuron{n}.con(c).weight);
                    else
                        newstr{s} = sprintf('%s,0)\n', newstr{s});
                        display('Caution: NetCon Weight initialized with default (0) !')
                    end
                    newstr{s} = sprintf('%sio = conList.append(con)',newstr{s});  %append con to conList
                    if nodeflag
                        newstr{s} = sprintf('%s}\n',newstr{s});
                    else
                        newstr{s} = sprintf('%s\n',newstr{s});
                    end
                end
                fprintf(ofile,strjoin(newstr));  % new connection
            end
            fprintf(ofile, 'objref con\n');
        end
        fprintf(ofile,'\n\n');
        fclose(ofile);
    elseif exist(fullfile(exchfolder,thisfolder,'init_con.hoc'),'file')
        delete(fullfile(exchfolder,thisfolder,'init_con.hoc'));
    end
    
    
    %% write init_rec.hoc
    
    x = getref(n,neuron,'record');
    x2 = getref(n,neuron,'APCount');
    x3 = getref(n,neuron,'pp');
    if (~isnan(x) || ~isnan(x2)) && (x==n || x2==n || x2~=x)     %rewrite only if one or both of record/APCount are not taken from previous sim or if both are taken from another sim but from different ones (not possible because both are in one hoc)
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_rec.hoc') ,'wt');   %open record hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define recording sites *****\n');
        if ~isnan(x)
            count = 0;  % counter for recording vector List
            countt = 0; % counter for time vector List
            for tt = 1:numel(thesetrees{n})
                t = thesetrees{n}(tt);
                if numel(neuron{x}.record) >= t && ~isempty(neuron{x}.record{t})  % if a recording site was defined for  this tree
                    recfields = fieldnames(neuron{x}.record{t});
                    if isfield(tree{thesetrees{n}(tt)},'artificial')
                        if numel(recfields) > 1 && strcmp(recfields,'record')
                            neuron{x}.record{t} = struct(tree{thesetrees{n}(tt)}.artificial,neuron{x}.record{t}); % put the structure in field named as the artificial neuron
                        end
                    end
                    recfields = fieldnames(neuron{x}.record{t});
                    
                    for f1 = 1:numel(recfields)
                        if isfield(tree{thesetrees{n}(tt)},'artificial')
                            rectype = 'artificial';
                        elseif strcmp(recfields{f1},'cell') %   size(neuron{x}.record{t},2)<3 || isempty(neuron{x}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                            rectype = 'cell';
                        else
                            rectype = 'pp';
                        end
                        
                        for r = 1:numel(neuron{x}.record{t}.(recfields{f1})) %.record)  % go through all variables to be recorded
                            
                            if strcmp(rectype,'cell')
                                if isfield(tree{thesetrees{n}(tt)},'R') && isfield(tree{thesetrees{n}(tt)},'rnames')
                                    Rs = tree{thesetrees{n}(tt)}.rnames(tree{thesetrees{n}(tt)}.R(neuron{x}.record{t}.cell(r).node));       % all region names of trees nodes
                                    strs = regexp(neuron{x}.record{t}.cell(r).record,'_','split');            % split record string to get mechanism name
                                    if numel(strs)>1   %any(strcmp(strs{1},{'v','i'}))             % check if record variable is variable of a mechanism or maybe global
                                        ignorethese = false(1,numel(neuron{x}.record{t}.cell(r).node));
                                        uRs = unique(Rs);
                                        str = '';
                                        for u = 1:numel(uRs)                                            % go through regions to be recorded
                                            if (~isfield(neuron{n}.mech{t},uRs{u}) || isfield(neuron{n}.mech{t},uRs{u}) && ~isfield(neuron{n}.mech{t}.(uRs{u}),strs{end})) &&  (~isfield(neuron{n}.mech{t},'all') || isfield(neuron{n}.mech{t},'all') && ~isfield(neuron{n}.mech{t}.all,strs{end}))             % check if this region also has the mechanism to be recorded
                                                ignorethese = ignorethese | strcmp(uRs{u},Rs);           % if not ignore these region for recording
                                                str = strcat(str,uRs{u},'/');
                                            end
                                        end
                                        if ~isempty(str)
                                            neuron{x}.record{t}.cell(r).node(ignorethese) = [];                     % delete the recording nodes which should be ignored
                                            warning(sprintf('Region(s) "%s" of tree %d do not contain mechanism "%s" for recording. Recording in this region is ignored',str(1:end-1),t,strs{end}))
                                        end
                                    end
                                end
                            end
                            
                            if ~any(strcmp(rectype,'artificial'))
                                inode = zeros(numel(neuron{x}.record{t}.(recfields{f1})(r).node),1);
                                for in = 1:numel(neuron{x}.record{t}.(recfields{f1})(r).node)
                                    inode(in) = find(minterf{thesetrees{n}(tt)}(:,1) == neuron{x}.record{t}.(recfields{f1})(r).node(in),1,'first');    %find the index of the node in minterf
                                end
                                [realrecs,~,ic] = unique(minterf{thesetrees{n}(tt)}(inode,[2,4]),'rows');
                                % put warning here !
                            end
                            
                            switch rectype
                                case 'cell'
                                    for in = 1:size(realrecs,1)
                                        fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('rec.label("%s at location %06.4f of section %d of cell %d")\n', neuron{x}.record{t}.cell(r).record , realrecs(in,2), realrecs(in,1) ,tt-1) ); % label the vector for plotting
                                        if params.cvode
                                            fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&%s(%f),rec,rect)}\n',tt-1,realrecs(in,1), neuron{x}.record{t}.cell(r).record, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{x}.record
                                        else
                                            fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),tvec)\n',tt-1,realrecs(in,1), neuron{x}.record{t}.cell(r).record, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{x}.record
                                        end
                                        
                                        fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                        if params.cvode
                                            fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                            neuron{x}.record{t}.cell(r).idt(in) = countt;   % reference to find recording in recList
                                            countt = countt +1;
                                        end
                                        neuron{x}.record{t}.cell(r).id(in) = count;  % reference to find recording in recList
                                        count = count +1;
                                    end
                                    neuron{x}.record{t}.cell(r).rrecs = realrecs; % gives the section and segment to the recordings
                                    neuron{x}.record{t}.cell(r).irrecs = ic; % gives the the index to realrecs for each node
                                case 'pp'
                                    delin = [];
                                    %                                 [~,iid] = intersect(neuron{x3}.pp{t}.(recfields{f1}).node,neuron{x}.record{t}.(recfields{f1}).node); % get reference to the node location of the PPs that should be connected
                                    for in =  1:size(realrecs,1)%numel(neuron{x}.record{t}{r,1})  % CAUTION might be wrong
                                        ind = find(neuron{x3}.pp{t}.(recfields{f1}).node == neuron{x}.record{t}.(recfields{f1})(r).node(in));
                                        if ~isempty(ind)
                                            fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                            fprintf(ofile,sprintf('rec.label("%s of %s Point Process at location %06.4f of section %d of cell %d")\n', neuron{x}.record{t}.(recfields{f1})(r).record , recfields{f1} , fliplr(realrecs(in,:)) ,tt-1) ); % label the vector for plotting
                                            if params.cvode
                                                fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                                fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&ppList.o(%d).%s,rec,rect)}\n',tt-1, realrecs(in,1),neuron{x3}.pp{t}.(recfields{f1})(r).id(ind), neuron{x}.record{t}.(recfields{f1})(r).record ) ); % record the parameter x at site y as specified in neuron{x}.record
                                            else
                                                fprintf(ofile,sprintf('io = rec.record(&ppList.o(%d).%s,tvec)\n',neuron{x3}.pp{t}.(recfields{f1})(r).id(ind), neuron{x}.record{t}.(recfields{f1})(r).record ) ); % record the parameter x at site y as specified in neuron{x}.record
                                            end
                                            fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                            if params.cvode
                                                fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                                neuron{x}.record{t}.(recfields{f1})(r).idt(in) = countt;   % reference to find recording in recList
                                                countt = countt +1;
                                            end
                                            neuron{x}.record{t}.(recfields{f1})(r).id(in) = count;   % reference to find recording in recList
                                            count = count +1;
                                        else
                                            delin = cat(1,delin,in);
                                            neuron{x}.record{t}.(recfields{f1})(r).id(in) = NaN;
                                            display(sprintf('Node %d of cell %d does not comprise the PP "%s". Recording is ignored.',neuron{x}.record{t}.(recfields{f1})(r).node(in),t,recfields{f1}))
                                            ic(ic == in) = [];  % if node does not correspond to some specified pp at that place, delete it
                                            ic(ic >= in) = ic(ic >= in) - 1;
                                        end
                                    end
                                    if ~isempty(delin)
                                        realrecs(delin,:) = [];  % if node does not correspond to some specified pp at that place, delete it
                                    end
                                    neuron{x}.record{t}.(recfields{f1})(r).rrecs = realrecs;
                                    neuron{x}.record{t}.(recfields{f1})(r).irrecs = ic; % gives the the index to realrecs for each node
                                case 'artificial'
                                    fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                    fprintf(ofile,sprintf('rec.label("%s of artificial cell %s (cell #%d)")\n', neuron{x}.record{t}.cell(r).record , tree{thesetrees{n}(tt)}.artificial, tt-1) ); % label the vector for plotting
                                    if strcmpi(neuron{x}.record{t}.cell(r).record,'on')
                                        fprintf(ofile,sprintf('nilcon = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',tt-1,0.5) );    % for art. cells, make netcon with threshold 0.5
                                        fprintf(ofile,sprintf('io = nilcon.record(rec)\n'));
                                        fprintf(ofile,'io = nilconList.append(nilcon)\n\n' );  %append recording vector to recList
                                        
                                    else
                                        if params.cvode
                                            fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                            fprintf(ofile,sprintf('io = cvode.record(&cellList.o(%d).cell.%s,rec,rect)\n',tt-1, neuron{x}.record{t}.cell(r).record ) );  % record the parameter x of artificial cell tt-1
                                        else
                                            fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).cell.%s,tvec)\n', tt-1, neuron{x}.record{t}.cell.record{r} ) ); % record the parameter x of artificial cell tt-1
                                        end
                                        if params.cvode
                                            fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                            neuron{x}.record{t}.cell(r).idt = countt;   % reference to find recording in recList
                                            countt = countt +1;
                                        end
                                    end
                                    fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                    neuron{x}.record{t}.cell(r).id = count;   % reference to find recording in recList
                                    count = count +1;
                            end
                        end
                        fprintf(ofile,'\n');
                    end
                end
            end
            fprintf(ofile, 'objref rec\n');
            fprintf(ofile, 'objref rect\n');
        end
        %
        if ~isnan(x2)     %rewrite only if record def is not taken from previous sim
            %!%! there might be problems with cvode (not adjusted yet)
            fprintf(ofile,'\n\n');
            fprintf(ofile,'// ***** Define APCount sites *****\n');
            %             if getref(n,neuron,'APCount') == n
            for tt = 1:numel(thesetrees{n})
                t = thesetrees{n}(tt);
                if numel(neuron{x2}.APCount) >= t && ~isempty(neuron{x2}.APCount{t})   % if a recording site was defined for  this tree
                    for r = 1: size(neuron{x2}.APCount{t},1)
                        if ~isfield(tree{thesetrees{n}(tt)},'artificial')
                            inode = find(minterf{thesetrees{n}(tt)}(:,1) == neuron{x2}.APCount{t}(r,1),1,'first');    %find the index of the node in minterf
                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',tt-1,minterf{thesetrees{n}(tt)}(inode,2) ) );    % corresponding section of node
                            fprintf(ofile,sprintf('{APC = new APCount(%f)\n',minterf{thesetrees{n}(tt)}(inode,3) ) );    % make APCCount at position x
                            fprintf(ofile,sprintf('APC.thresh = %f\n',neuron{x2}.APCount{t}(r,2) ) ); % set threshold of APCount [mV]
                        else
                            fprintf(ofile,sprintf('APC = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',tt-1,neuron{x2}.APCount{t}(r,2) ) );    % for art. cells, make netcon with threshold
                        end
                        fprintf(ofile,'APCrec = new Vector()\n');
                        fprintf(ofile,'io = APCrecList.append(APCrec)\n');
                        fprintf(ofile,'io = APC.record(APCrecList.o(APCrecList.count()-1))\n');
                        
                        if ~isfield(tree{thesetrees{n}(tt)},'artificial')
                            fprintf(ofile,'io = APCList.append(APC)}\n\n' );  %append recording vector to recList
                        else
                            fprintf(ofile,'io = APCList.append(APC)\n\n' );  %append recording vector to recList
                        end
                    end
                    fprintf(ofile,'\n');
                end
            end
            fprintf(ofile, 'objref APC\n');
            fprintf(ofile, 'objref APCrec\n');
            %             end
        end
        fclose(ofile);
    elseif exist(fullfile(exchfolder,thisfolder,'init_rec.hoc'),'file')
        delete(fullfile(exchfolder,thisfolder,'init_rec.hoc'));
    end
    
    
    %% write init_play.hoc
    
    x = getref(n,neuron,'play');
    if (~isnan(x)) && x==n      %rewrite only if one or both of play/APCount are not taken from previous sim or if both are taken from another sim but from different ones (not possible because both are in one hoc)
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_play.hoc') ,'wt');   %open play hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define play sites *****\n');
        %     end
        count = 0;  % counter for playing vector List
        countt = 0; % counter for time vector List
        for tt = 1:numel(thesetrees{n})
            t = thesetrees{n}(tt);
            if numel(neuron{x}.play) >= t && ~isempty(neuron{x}.play{t})  % if a playing site was defined for  this tree
                
                if isfield(tree{thesetrees{n}(tt)},'artificial')
                    playfields = fieldnames(neuron{x}.play{t});
                    if numel(playfields) > 1 && strcmp(playfields,'play')
                        neuron{x}.play{t} = struct(tree{thesetrees{n}(tt)}.artificial,neuron{x}.play{t}); % put the structure in field named as the artificial neuron
                    end
                end
                playfields = fieldnames(neuron{x}.play{t});
                
                for f1 = 1:numel(playfields)
                    if isfield(tree{thesetrees{n}(tt)},'artificial')
                        playtype = 'artificial';
                    elseif strcmp(playfields{f1},'cell') %   size(neuron{x}.play{t},2)<3 || isempty(neuron{x}.play{t}{r,3})       % check if playing should be a parameter in a section, or a point process (divided in pps and electrodes)
                        playtype = 'cell';
                    else
                        playtype = 'pp';
                    end
                    
                    for r = 1:numel(neuron{x}.play{t}.(playfields{f1})) %.play)  % go through all variables to be played
                        if isfield(neuron{n}.play{t}.(playfields{f1})(r),'cont')
                            cont = neuron{n}.play{t}.(playfields{f1})(r).cont;
                        else
                            cont = 0;
                        end
                        if strcmp(playtype,'cell')
                            if isfield(tree{thesetrees{n}(tt)},'R') && isfield(tree{thesetrees{n}(tt)},'rnames')
                                Rs = tree{thesetrees{n}(tt)}.rnames(tree{thesetrees{n}(tt)}.R(neuron{x}.play{t}.cell(r).node));       % all region names of trees nodes
                                strs = regexp(neuron{x}.play{t}.cell(r).play,'_','split');            % split play string to get mechanism name
                                if numel(strs)>1   %any(strcmp(strs{1},{'v','i'}))             % check if play variable is variable of a mechanism or maybe global
                                    ignorethese = false(1,numel(neuron{x}.play{t}.cell(r).node));
                                    uRs = unique(Rs);
                                    str = '';
                                    for u = 1:numel(uRs)                                            % go through regions to be played
                                        if (~isfield(neuron{n}.mech{t},uRs{u}) || isfield(neuron{n}.mech{t},uRs{u}) && ~isfield(neuron{n}.mech{t}.(uRs{u}),strs{end})) &&  (~isfield(neuron{n}.mech{t},'all') || isfield(neuron{n}.mech{t},'all') && ~isfield(neuron{n}.mech{t}.all,strs{end}))             % check if this region also has the mechanism to be played
                                            ignorethese = ignorethese | strcmp(uRs{u},Rs);           % if not ignore these region for playing
                                            str = strcat(str,uRs{u},'/');
                                        end
                                    end
                                    if ~isempty(str)
                                        neuron{x}.play{t}.cell(r).node(ignorethese) = [];                     % delete the playing nodes which should be ignored
                                        warning(sprintf('Region(s) "%s" of tree %d do not contain mechanism "%s" for playing. Playing in this region is ignored',str(1:end-1),t,strs{end}))
                                    end
                                end
                            end
                        end
                        
                        if ~any(strcmp(playtype,'artificial'))
                            inode = zeros(numel(neuron{x}.play{t}.(playfields{f1})(r).node),1);
                            for in = 1:numel(neuron{x}.play{t}.(playfields{f1})(r).node)
                                inode(in) = find(minterf{thesetrees{n}(tt)}(:,1) == neuron{x}.play{t}.(playfields{f1})(r).node(in),1,'first');    %find the index of the node in minterf
                            end
                            [realplays,~,ic] = unique(minterf{thesetrees{n}(tt)}(inode,[2,4]),'rows');
                            % put warning here !
                        end
                        
                        switch playtype
                            case 'cell'
                                for in = 1:size(realplays,1)
                                    
                                    fprintf(ofile,sprintf('playt = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).times) ));    % create new playing time vector
                                    %a file needs to be created to temporally save the vector so
                                    %NEURON can read it in. otherwise it would be necessary to
                                    %print the whole vector into the hoc file. alternatively i
                                    %could give a file name where the vector lies so it is not
                                    %written each time cn is called...
                                    f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_at_%d_cell_%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1)),'w');
                                    fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).times(1:end-1));
                                    fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).times(end));
                                    fclose(f);
                                    fprintf(ofile,'f = new File()\n');
                                    fprintf(ofile,sprintf('f.ropen("plt_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in),tt-1));  %vector file is opened
                                    fprintf(ofile,'playt.scanf(f)\n');    % file is read into time vector
                                    fprintf(ofile,'io = f.close()\n');     %file is closed
                                    fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append playing time vector to playtList
                                    
                                    fprintf(ofile,sprintf('play = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).value) ) );    % create new playing vector
                                    f = fopen(fullfile(exchfolder,thisfolder,sprintf('pl_%s_at_%d_cell_%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1)),'w');
                                    fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).value(1:end-1));
                                    fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).value(end));
                                    fclose(f);
                                    fprintf(ofile,'f = new File()\n');
                                    fprintf(ofile,sprintf('f.ropen("pl_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1));  %vector file is opened
                                    fprintf(ofile,'play.scanf(f)\n');     % file is read into play vector
                                    fprintf(ofile,'io = f.close()\n');   %file is closed
                                    fprintf(ofile,sprintf('play.label("playing %s at node %d of cell %d")\n', neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in) ,tt-1) ); % label the vector for plotting
                                    fprintf(ofile,sprintf('play.play(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),playtList.o(playtList.count()-1),%d)\n',tt-1,realplays(in,1), neuron{n}.play{t}.(playfields{f1})(r).play, realplays(in,2), cont ) ); % play the parameter x at site y as specified in neuron{n}.play
                                    fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
                                    
                                end
                                %                                 neuron{x}.play{t}.cell(r).rplays = realplays; % gives the section and segment to the playings
                                %                                 neuron{x}.play{t}.cell(r).irplays = ic; % gives the the index to realplays for each node
                            case 'pp'
                                %                                 delin = [];
                                for in =  1:numel(inode)%size(realplays,1)%numel(neuron{x}.play{t}{r,1})  % CAUTION might be wrong
                                    ind = find(neuron{x3}.pp{t}.(playfields{f1}).node == neuron{x}.play{t}.(playfields{f1})(r).node(in));  % check if node really has the corresponding PP
                                    if ~isempty(ind)
                                        fprintf(ofile,sprintf('playt = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).times) ));    % create new playing time vector
                                        %a file needs to be created to temporally save the vector so
                                        %NEURON can read it in. otherwise it would be necessary to
                                        %print the whole vector into the hoc file. alternatively i
                                        %could give a file name where the vector lies so it is not
                                        %written each time cn is called...
                                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_%s_at_%d_cell_%d.dat', playfields{f1}, neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1)),'w');
                                        if ~any(size(neuron{n}.play{t}.(playfields{f1})(r).times)==1) % it's a matrix
                                            if size(neuron{n}.play{t}.(playfields{f1})(r).times,1) == numel(inode)
                                                fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).times(in,1:end-1));
                                                fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).times(in,end));
                                            else
                                                for nn = 1:numel(neuron)
                                                    delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                end
                                                error('Times vector of play feature has wrong size')
                                            end
                                        else
                                            fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).times(1:end-1));
                                            fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).times(end));
                                        end
                                        fclose(f);
                                        fprintf(ofile,'f = new File()\n');
                                        fprintf(ofile,sprintf('f.ropen("plt_%s_%s_at_%d_cell_%d.dat")\n', playfields{f1}, neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in),tt-1));  %vector file is opened
                                        fprintf(ofile,'playt.scanf(f)\n');    % file is read into time vector
                                        fprintf(ofile,'io = f.close()\n');     %file is closed
                                        fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append playing time vector to playtList
                                        
                                        fprintf(ofile,sprintf('play = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).value) ) );    % create new playing vector
                                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('pl_%s_%s_at_%d_cell_%d.dat', playfields{f1}, neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1)),'w');
                                        fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).value(1:end-1));
                                        fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).value(end));
                                        fclose(f);
                                        fprintf(ofile,'f = new File()\n');
                                        fprintf(ofile,sprintf('f.ropen("pl_%s_%s_at_%d_cell_%d.dat")\n', playfields{f1}, neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,tt-1));  %vector file is opened
                                        fprintf(ofile,'play.scanf(f)\n');     % file is read into play vector
                                        fprintf(ofile,'io = f.close()\n');   %file is closed
                                        fprintf(ofile,sprintf('play.label("playing %s %s at node %d of cell %d")\n', playfields{f1}, neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in) ,tt-1) ); % label the vector for plotting
                                        fprintf(ofile,sprintf('io = play.play(&ppList.o(%d).%s,playt)\n',neuron{x3}.pp{t}.(playfields{f1})(r).id(ind), neuron{x}.play{t}.(playfields{f1})(r).play ) ); % play the parameter x at site y as specified in neuron{x}.play
                                        fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
                                        
                                        neuron{x}.play{t}.(playfields{f1})(r).id(in) = count;   % reference to find playing in playList
                                        count = count +1;
                                    else
                                        %                                         delin = cat(1,delin,in);
                                        neuron{x}.play{t}.(playfields{f1})(r).id(in) = NaN;   % reference to find playing in playList
                                        display(sprintf('Node %d of cell %d does not comprise the PP "%s". Playing is ignored.',neuron{x}.play{t}.(playfields{f1})(r).node(in),t,playfields{f1}))
                                        %                                         ic(ic == in) = [];  % if node does not correspond to some specified pp at that place, delete it
                                        %                                         ic(ic >= in) = ic(ic >= in) - 1;
                                    end
                                    %                                     if ~isempty(delin)
                                    %                                         realplays(delin,:) = [];  % if node does not correspond to some specified pp at that place, delete it
                                    %                                     end
                                end
                                %                                 neuron{x}.play{t}.(playfields{f1})(r).rplays = realplays;
                                %                                 neuron{x}.play{t}.(playfields{f1})(r).irplays = ic; % gives the the index to realplays for each node
                            case 'artificial'
                                fprintf(ofile,sprintf('playt = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).times) ));    % create new playing time vector
                                %a file needs to be created to temporally save the vector so
                                %NEURON can read it in. otherwise it would be necessary to
                                %print the whole vector into the hoc file. alternatively i
                                %could give a file name where the vector lies so it is not
                                %written each time cn is called...
                                if strcmp(tree{thesetrees{n}(tt)}.artificial,'VecStim') && any(neuron{n}.play{t}.(playfields{f1})(r).times < 0)
                                    neuron{n}.play{t}.(playfields{f1})(r).times(neuron{n}.play{t}.(playfields{f1})(r).times < 0) = [];
                                    display('Warning, VecStim should not receive negative play times. These are deleted now')
                                end
                                f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_of_art_%s_cell%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play, playfields{f1}, tt-1)),'w');
                                fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).times(1:end-1));
                                fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).times(end));
                                fclose(f);
                                fprintf(ofile,'f = new File()\n');
                                fprintf(ofile,sprintf('f.ropen("plt_%s_of_art_%s_cell%d.dat")\n',  neuron{n}.play{t}.(playfields{f1})(r).play, playfields{f1}, tt-1));  %vector file is opened
                                fprintf(ofile,'playt.scanf(f)\n');    % file is read into time vector
                                fprintf(ofile,'io = f.close()\n');     %file is closed
                                fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append playing time vector to playtList
                                %                                     end
                                if ~strcmp(neuron{n}.play{t}.(playfields{f1})(r).play,'spike') && ~strcmp(playfields{f1},'VecStim')
                                    fprintf(ofile,sprintf('play = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).value) ) );    % create new playing vector
                                    f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_of_art_%s_cell%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play, playfields{f1}, tt-1)),'w');
                                    fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).value(1:end-1));
                                    fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).value(end));
                                    fclose(f);
                                    fprintf(ofile,'f = new File()\n');
                                    fprintf(ofile,sprintf('f.ropen("plt_%s_of_art_%s_cell%d.dat")\n', neuron{n}.play{t}.(playfields{f1})(r).play, playfields{f1}, tt-1));  %vector file is opened
                                    fprintf(ofile,'play.scanf(f)\n');     % file is read into play vector
                                    fprintf(ofile,'io = f.close()\n');   %file is closed
                                    fprintf(ofile,sprintf('play.label("%s of artificial cell %s (cell #%d)")\n', neuron{x}.play{t}.cell(r).play , tree{thesetrees{n}(tt)}.artificial, tt-1) ); % label the vector for plotting                                    fprintf(ofile,sprintf('io = play.play(&ppList.o(%d).%s,playt)\n',neuron{x3}.pp{t}.(playfields{f1})(r).id(ind), neuron{x}.play{t}.(playfields{f1})(r).play ) ); % play the parameter x at site y as specified in neuron{x}.play
                                    warning('not tested yet play')
                                    fprintf(ofile,sprintf('io = cellList.o(%d).cell.play(&,playt)\n',tt-1, neuron{x}.play{t}.(playfields{f1})(r).play ) ); % play the parameter x at site y as specified in neuron{x}.play
                                    fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
                                else
                                    fprintf(ofile,sprintf('io = cellList.o(%d).cell.play(playt)\n',tt-1 ) ); % the vecstim has its own play class
                                end
                                
                                neuron{x}.play{t}.cell(r).id = count;   % reference to find playing in playList
                                count = count +1;
                                
                        end
                    end
                    fprintf(ofile,'\n');
                end
            end
        end
        fprintf(ofile, 'objref play\n');
        fprintf(ofile, 'objref playt\n');
        %
        fclose(ofile);
    elseif exist(fullfile(exchfolder,thisfolder,'init_play.hoc'),'file')
        delete(fullfile(exchfolder,thisfolder,'init_play.hoc'));
    end
    
    %% write save_rec.hoc
    
    ofile = fopen(fullfile(exchfolder,thisfolder,'save_rec.hoc') ,'wt');   %open record hoc file in write modus
    fprintf(ofile,'// * Write Recordings to Files *\n');
    %     end
    x = getref(n,neuron,'record');
    x2 = getref(n,neuron,'APCount');
    x3 = getref(n,neuron,'pp');
    if ~isnan(x)
        out{n}.record = cell(1,numel(thesetrees{n}));   % initialize output of cn
        
        for tt = 1:numel(thesetrees{n})
            t = thesetrees{n}(tt);
            if numel(neuron{x}.record) >= t && ~isempty(neuron{x}.record{t})
                %                 for r = 1: size(neuron{x}.record{t},1)
                recfields = fieldnames(neuron{x}.record{t});
                
                for f1 = 1:numel(recfields)
                    if isfield(tree{thesetrees{n}(tt)},'artificial')
                        rectype = 'artificial';
                    elseif strcmp(recfields{f1},'cell') %   size(neuron{x}.record{t},2)<3 || isempty(neuron{x}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                        rectype = 'cell';
                    else
                        rectype = 'pp';
                    end
                    
                    for r = 1:numel(neuron{x}.record{t}.(recfields{f1})) %.record)  % go through all variables to be recorded
                        
                        
                        switch rectype
                            case 'cell'
                                for in = 1:size(neuron{x}.record{t}.cell(r).rrecs,1)
                                    fname = sprintf('cell%d_sec%d_loc%06.4f_%s', tt-1, neuron{x}.record{t}.cell(r).rrecs(in,:), neuron{x}.record{t}.(recfields{f1})(r).record  );
                                    fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                    fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname) )  );  % open file for this vector with write perm.
                                    fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).id(in) ) );    % print the data of the vector into the file
                                    fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                    if params.cvode
                                        fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                        fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname) )  );  % open file for this vector with write perm.
                                        fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).idt(in) ) );    % print the data of the vector into the file
                                        fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                    end
                                    noutfiles = noutfiles +1;
                                    %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',tt-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
                                    %                                     readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{x}.record{t}.(recfields{f1}).record{r} , neuron{x}.record{t}.(recfields{f1}).rrecs(in,1), neuron{x}.record{t}.(recfields{f1}).rrecs(in,2) };
                                    readfiles{noutfiles} = {sprintf('%s.dat',fname), n, t , 'cell', neuron{x}.record{t}.(recfields{f1})(r).record , neuron{x}.record{t}.(recfields{f1})(r).node(neuron{x}.record{t}.(recfields{f1})(r).irrecs == in) }; %neuron{x}.record{t}.(recfields{f1})(r).node(in) };
                                end
                            case 'pp'
                                for in = 1:size(neuron{x}.record{t}.(recfields{f1})(r).rrecs,1)
                                    if ~isnan(neuron{x}.record{t}.(recfields{f1})(r).id(in))  % if recording has not been deleted because the PP did not exist at that place
                                        fname = sprintf('%s_cell%d_sec%d_loc%06.4f_%s',recfields{f1},tt-1, neuron{x}.record{t}.(recfields{f1})(r).rrecs(in,:), neuron{x}.record{t}.(recfields{f1})(r).record  );
                                        
                                        fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                        fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname))  );  % open file for this vector with write perm.
                                        fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).id(in) ) );    % print the data of the vector into the file
                                        fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                        if params.cvode
                                            fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                            fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname))  );  % open file for this vector with write perm.
                                            fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).idt(in) ) );    % print the data of the vector into the file
                                            fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                        end
                                        noutfiles = noutfiles +1;
                                        %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',tt-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
                                        %                                     readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{x}.record{t}.(recfields{f1}).record{r} };
                                        readfiles{noutfiles} = {sprintf('%s.dat',fname), n, t , recfields{f1}, neuron{x}.record{t}.(recfields{f1})(r).record , neuron{x}.record{t}.(recfields{f1})(r).node(neuron{x}.record{t}.(recfields{f1})(r).irrecs == in)};%neuron{x}.record{t}.(recfields{f1})(r).node(in) };
                                    end
                                end
                            case 'artificial'
                                fname = sprintf('cell%d_%s',tt-1, neuron{x}.record{t}.(recfields{f1})(r).record );
                                fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname))  );  % open file for this vector with write perm.
                                fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).id ) );    % print the data of the vector into the file
                                fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                if params.cvode
                                    fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                    fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname))  );  % open file for this vector with write perm.
                                    fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', neuron{x}.record{t}.(recfields{f1})(r).idt ) );    % print the data of the vector into the file
                                    fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                end
                                noutfiles = noutfiles +1;
                                %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',tt-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
                                %                                 readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{x}.record{t}.(recfields{f1}).record{r} };
                                readfiles{noutfiles} = {sprintf('%s.dat',fname), n, t , 'cell', neuron{x}.record{t}.(recfields{f1})(r).record , 1 };
                        end
                    end
                end
            end
        end
        fprintf(ofile,'\n');
    end
    fprintf(ofile,'// * Write APCounts to Files *\n');
    
    x = getref(n,neuron,'APCount');
    if ~isnan(x)
        c=0;
        for tt = 1:numel(thesetrees{n})
            t = thesetrees{n}(tt);
            if numel(neuron{x}.APCount) >= t && ~isempty(neuron{x}.APCount{t})     % if a recording site was defined for  this tree
                for r = 1: size(neuron{x}.APCount{t},1)
                    fname = sprintf('cell%d_node%d_APCtimes.dat',tt-1,neuron{x}.APCount{t}(r,1) );
                    fprintf(ofile,'f = new File()\n');      %create a new filehandle
                    fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,fname) );  % open file for this vector with write perm.
                    fprintf(ofile,sprintf('io = APCrecList.o(%d).printf(f, "%%%%-20.10g")\n', c ) );    % print the data of the vector into the file
                    fprintf(ofile,'io = f.close()\n');   %close the filehandle
                    
                    c= c+1;
                    noutfiles = noutfiles +1;
                    %                     readfiles{noutfiles} = {n, fname , 'APCtimes' , t , neuron{x}.APCount{t}(r,1) };
                    readfiles{noutfiles} = {fname, n, t , 'APCtimes', 'times' , neuron{x}.APCount{t}(r,1)};
                end
                fprintf(ofile,'\n');
            end
        end
    end
    fclose(ofile);
    
    
    if ~isempty(strfind(options,'-cl')) %transfer files to server
        filenames = {interf_file,'init_cells.hoc','init_mech.hoc','init_pp.hoc','init_con.hoc','init_rec.hoc','save_rec.hoc','init_play.hoc'}; %'init_pas.hoc','init_stim.hoc'
        m = 1;
        localfilename{m} = fullfile(exchfolder,thisfolder,filenames{1});
        remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{1});
        m = m + 1;
        if usestreesof(n) == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{2});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{2});
            m = m + 1;
        end
        if  getref(n,neuron,'mech') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{3});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{3});
            m = m + 1;
        end
        if getref(n,neuron,'pp') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{4});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{4});
            m = m + 1;
        end
        if getref(n,neuron,'con') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{5});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{5});
            m = m + 1;
        end
        if getref(n,neuron,'record') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{6});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{6});
            m = m + 1;
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{7});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{7});
            m = m + 1;
        end
        if getref(n,neuron,'play') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{8});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{8});
            m = m + 1;
        end
            %create job
            ofile = fopen(fullfile(exchfolder,thisfolder,'start_nrn.pbs') ,'wt');
            fprintf(ofile,'#!/bin/bash\n');
            fprintf(ofile,'# set job variables\n');
            fprintf(ofile,'#$ -j n\n');
            fprintf(ofile,'#$ -S /bin/bash\n');
            fprintf(ofile,'#$ -pe openmp 1\n');
            fprintf(ofile,sprintf('#$ -l h_rt=%02d:%02d:%02d\n',params.server.walltime));
            fprintf(ofile,sprintf('#$ -l h_rt=%02d:%02d:%02d\n',params.server.softwalltime));
            fprintf(ofile,sprintf('#$ -l h_vmem=%02dG\n',params.server.memory));
            fprintf(ofile,sprintf('#$ -N nrn_%s\n',thisfolder));  % name of job
%             fprintf(ofile,sprintf('#$ -o %s/%s/%s',nrn_exchfolder,thisfolder,'std.out'));
%             fprintf(ofile,sprintf('#$ -e %s/%s/%s',nrn_exchfolder,thisfolder,'err.out'));
            fprintf(ofile,'# load needed modules \n');
            fprintf(ofile,sprintf('module load %s \n ',params.server.neuron{1})); % load first found neuron version
%             fprintf(ofile,'module load openmpi/gcc/64/1.3.3\n');
%             fprintf(ofile,'# change to path with your executable\n');
                
            fprintf(ofile,'# start program\n');
            fprintf(ofile,sprintf('$NRNHOME/x86_64/nrniv -nobanner -nogui "%s/%s/%s" > "%s/%s/NeuronLogFile.txt" 2> "%s/%s/ErrorLogFile.txt" \n',nrn_exchfolder,thisfolder,interf_file,nrn_exchfolder,thisfolder,nrn_exchfolder,thisfolder));

            fclose(ofile);
            localfilename{m} = fullfile(exchfolder,thisfolder,'start_nrn.pbs');
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,'start_nrn.pbs');
%         params.server.connect = sftpfrommatlab(params.server.connect,localfilename,remotefilename);
        sftpfrommatlab(params.server.user,params.server.host,params.server.pw,localfilename,remotefilename);
    end
    
    if strfind(options,'-d')
        tim = toc(tim);
        fprintf(sprintf('Sim %d: HOC writing time: %g min %.2f sec\n',n,floor(tim/60),rem(tim,60)))
    end
    
end

%% Execute NEURON
if ~isempty(strfind(options,'-cl'))
    num_cores = 500;  % use evt qstat -f?
else
    num_cores = feature('numCores');
end
if isfield(params,'numCores')
    if num_cores < params.numCores
        warning(sprintf('%d cores have been assigned to T2N, however only %d physical cores where detected. Defining more cores might slow down PC and simulations',params.numCores,num_cores))
    end
    num_cores = params.numCores;
    
end
    
simids = zeros(numel(neuron),1); % 0 = not run, 1 = running, 2 = finished, 3 = error
jobid = simids;
flag = logical(jobid);
for s = 1:num_cores
    r = find(simids==0,1,'first');
    if ~isempty(r)
        [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);
        simids(r) = 1;
    end
end

if noutfiles > 0 % if output is expected
    % wait for the NEURON process to be finished as indicated by appearance of
    % a file called 'readyflag' in the exchfolder; should be created in the last
    % line of the NEURON program
    if isempty(strfind(options,'-q'))
        display('waiting for NEURON to finish...')
    end
    if ~isempty(strfind(options,'-w'))
        w = waitbar(0,'Neuron Simulations are running, please wait');
    end
    timm = tic;
    while ~all(simids>1)
        if ~isempty(strfind(options,'-cl'))
                pause(0.8);
                [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('%s qstat -u %s',params.server.envstr,params.server.user));
                for s = find(simids==1)
%                     pause(0.2);
                    
                    if ~isempty(regexp(answer.StdOut,jobid(s),'ONCE'))      %job not finished yet
                        regexp(answer.StdOut,jobid(s))
                        answer = textscan(answer.StdOut,'%*[^QR]%s%*[^QR]');
                        if strcmpi(answer.StdOut,'R')
                            if ~flag(s)
                                display(sprintf('Simulation %d is still running on cluster',s))
                                if strfind(options,'-d')
                                    tim = toc(tim);
                                    display(sprintf('Cluster queue wait time: %g min %.2f sec',floor(tim/60),rem(tim,60)))
                                    tim = tic;  %reset timer to acount for queue wait time
                                    flag(s) = true;
                                end
                            end
                            %                     pause((params.server.walltime(1)*60+params.server.walltime(2))*60+params.server.walltime(3))
                        end
                    else   % job is finished
                        simids(s) = 2;              % mark that simulation as finished
                        r = find(simids==0,1,'first');  % find next not runned simid
                        [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                        simids(r) = 1;          % mark this as running
                    end
                    %             [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/readyflag',nrn_exchfolder));
                end
        else
            s = find(simids==1);
            for ss = 1:numel(s)
                r = exist(fullfile(exchfolder,sprintf('sim%d',s(ss)),'readyflag'),'file');  % becomes 1 (still running) if not existing, or 2 (finished)
                if r == 2
                    simids(s(ss)) = 2;              % mark that simulation as finished
                    r = find(simids==0,1,'first');  % find next not runned simid
                    if ~isempty(r)
                        [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                        simids(r) = 1;          % mark this as running
                    end
                elseif exist(fullfile(exchfolder,sprintf('sim%d',s(ss)),'ErrorLogFile.txt'),'file') == 2
                    finfo = dir(fullfile(exchfolder,sprintf('sim%d',s(ss)),'ErrorLogFile.txt'));
                    if finfo.bytes > 0      % because error file log is always built
                        f = fopen(fullfile(exchfolder,sprintf('sim%d',s(ss)),'ErrorLogFile.txt'));
                        txt = fscanf(f,'%c');
                        fclose(f);
                        errordlg(sprintf('There was an error in Simulation %d:\n******************************\n%s\n******************************\nDue to that t2n has no output to that Simulation.',s(ss),txt(1:min(numel(txt),2000))));
                        simids(s(ss)) = 3;
                        r = find(simids==0,1,'first');  % find next not runned simid
                        if ~isempty(r)
                            [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                            simids(r) = 1;          % mark this as running
                        end
                    end
                end
            end
            pause(0.1);
        end
        if ~isempty(strfind(options,'-w'))
            if ishandle(w)
                waitbar(sum(simids>1)/numel(simids),w)
            else
                answer = questdlg(sprintf('Waitbar was closed, t2n stopped continuing. Only finished data is returned. If accidently, retry.\nClose all NEURON instances?\n (Caution if several Matlab instances are running)'),'Close NEURON instances?','Close','Ignore','Ignore');
                if strcmp(answer,'Close')
                    system('taskkill /F /IM nrniv.exe');
                end
                simids(simids<2) = 4;
                fclose all;
            end
        end
    end
    
    if strfind(options,'-cl')
        s = find(simids==1);
        timm = tic;
        for ss = 1:numel(s)
            [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/%s/readyflag',nrn_exchfolder,sprintf('sim%d',s(ss))));
            if isempty(answer.StdOut)    % then there was an error during executing NEURON
                [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('cat %s/%s/ErrorLogFile.txt',nrn_exchfolder,sprintf('sim%d',s(ss))));
%                 bash-4.2$ cat wait.out | tail
                if ~isempty(answer.StdOut)
                    error('There was an error during NEURON simulation:\n%s\n',answer.StdOut)
                end
                r = find(simids==0,1,'first');  % find next not runned simid
                [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                simids(r) = 1;          % mark this as running
                simids(s(ss)) = 3;
                %                 return
            end
        end
    end
    
    if ~isempty(strfind(options,'-d'))
        tim = toc(timm);
        fprintf(sprintf('NEURON execute time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
    if isempty(strfind(options,'-q'))
        display('NEURON finished... loading data...')
    end
    if ~isempty(strfind(options,'-w')) && ishandle(w)
        close(w)
    end
    if strfind(options,'-d')
        tim = tic;
    end
    if strfind(options,'-cl')
        outputnames = cellfun(@(y) strcat(nrn_exchfolder,sprintf('/sim%d/',y{1}),y{2}),readfiles,'UniformOutput',0);  % extract filenames
        scptomatlab(params.server.connect,exchfolder,outputnames)
    end
    
    
    if ~params.cvode
        % load time vector from NEURON (necessary because of roundoff errors
        for s = 1:numel(simids)
            fn = fullfile(exchfolder,sprintf('sim%d',s),'tvec.dat');
            out{s}.t = load(fn,'-ascii');
        end
    end
    
    %% Receive files from Neuron
    if ~isempty(strfind(options,'-w'))
        w = waitbar(0,'Loading files, please wait');
    end
    for f = 1:noutfiles
        if simids(readfiles{f}{2}) == 2    % if there was no error during simulation
            fn = fullfile(exchfolder,sprintf('sim%d',readfiles{f}{2}),readfiles{f}{1});
            if numel(readfiles{f}{6}) > 1
                sprintf('Recording of %s in %s has redundant values since nodes are in same segment.\n',readfiles{f}{5},readfiles{f}{4})
            end
            switch readfiles{f}{4}
                case 'APCtimes'
                    out{readfiles{f}{2}}.APCtimes{readfiles{f}{3}}(readfiles{f}{6}) = repmat({load(fn,'-ascii')},numel(readfiles{f}{6}),1);
                otherwise
                    out{readfiles{f}{2}}.record{readfiles{f}{3}}.(readfiles{f}{4}).(readfiles{f}{5})(readfiles{f}{6}) = repmat({load(fn,'-ascii')},numel(readfiles{f}{6}),1);
                    if params.cvode
                        if params.use_local_dt  % if yes, dt was different for each cell, so there is more than one time vector
                            if numel(out{readfiles{f}{2}}.t) < readfiles{f}{3} || isempty(out{readfiles{f}{2}}.t{readfiles{f}{3}})
                                out{readfiles{f}{2}}.t{readfiles{f}{3}} = load(strcat(fn(1:end-4),'_tvec',fn(end-3:end)),'-ascii');    %loading of one time vector file per cell (sufficient)
                            end
                        elseif ~ isfield(out{readfiles{f}{2}},'t')       % if it has been loaded in a previous loop
                            out{readfiles{f}{2}}.t = load(strcat(fn(1:end-4),'_tvec',fn(end-3:end)),'-ascii');    %loading of one time vector file at all (sufficient)
                        end
                        out{readfiles{f}{2}}.t(find(diff(out{readfiles{f}{2}}.t,1) == 0) + 1) = out{readfiles{f}{2}}.t(find(diff(out{readfiles{f}{2}}.t,1) == 0) + 1) + 1e-10;  % add tiny time step to tvec to avoid problems with step functions
                    end
                    
            end
            delete(fn)  % delete dat file after loading
        elseif simids(readfiles{f}{2}) == 4  % t2n was aborted
            out{readfiles{f}{2}}.error = 2;
        else
            out{readfiles{f}{2}}.error = 1;
        end
        if ~isempty(strfind(options,'-w'))
            if ishandle(w)
                waitbar(f/noutfiles,w);
            else
                answer = questdlg(sprintf('Waitbar has been closed during data loading. If accidently, retry.\nClose all NEURON instances?\n (Caution if several Matlab instances are running)'),'Close NEURON instances?','Close','Ignore','Ignore');
                if strcmp(answer,'Close')
                    system('taskkill /F /IM nrniv.exe');
                end
                fclose all;
                for n = 1:numel(neuron)
                    delete(fullfile(exchfolder,sprintf('sim%d',n),'iamrunning'));   % delete the running mark
                end
                return
            end
            
        end
    end
    
    if ~params.cvode
        for n = 1:numel(neuron)
            out{n}.t = load(fullfile(exchfolder,sprintf('sim%d',n),'tvec.dat'));
        end
    end
    
    if isempty(strfind(options,'-q'))
        display('data sucessfully loaded')
    end
    if ~isempty(strfind(options,'-w'))
        close(w)
    end
    if strfind(options,'-d')
        tim = toc(tim);
        fprintf(sprintf('Data loading time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
    
end
if outoptions.nocell
    out = out{1};
end
for n = 1:numel(neuron)
    delete(fullfile(exchfolder,sprintf('sim%d',n),'iamrunning'));   % delete the running mark
end


end


function [jobid,tim] = exec_neuron(simid,exchfolder,nrn_exchfolder,interf_file,params,options)
%% Execute NEURON
tim = tic;

% switch params.openNeuron
%     case 1
%         opsign = '&';%sprintf(' > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt" &',exchfolder,simid,exchfolder,simid); % execute the program in foreground and hand control to NEURON
%     case 0
%         opsign = sprintf(' -c quit() > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt" \n exit',exchfolder,simid,exchfolder,simid);  % execute the program iconified
% end
% execute the file in neuron:
fname = regexprep(fullfile(exchfolder,sprintf('sim%d',simid),interf_file),'\\','/');

if ~isempty(strfind(options,'-cl'))
       [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('%s%sqsub -p 0 %s/%s/%s',params.server.envstr,params.server.qfold,nrn_exchfolder,sprintf('sim%d',simid),'start_nrn.pbs'));
%        [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('%s/%s/%s',nrn_exchfolder,sprintf('sim%d',simid),'start_nrn.job'));
        fprintf(sprintf('Answer server after submitting:\n%s\n%s\nExtracing Job Id and wait..\n',answer.StdOut,answer.StdErr))
else
    oldpwd = '';
    if ~strcmpi(pwd,params.path)
        oldpwd = pwd;
        cd(params.path);   % in order to have NEURON the path as its starting folder
    end
    if params.openNeuron
        %         dos([params.neuronpath ' -nobanner "' fname '" ' opsign]); %&,char(13),'exit&']); %nrniv statt neuron
        system(['start ' params.neuronpath ' -nobanner "' fname sprintf('" > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt"',exchfolder,simid,exchfolder,simid)]); %&,char(13),'exit&']); %nrniv statt neuron
    else
        system(['start /B ' params.neuronpath ' -nobanner "' fname sprintf('" -c quit() > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt"',exchfolder,simid,exchfolder,simid)]); %&,char(13),'exit&']); %nrniv statt neuron
    end
    %         system(['wmic process call create ''', params.neuronpath, ' -nobanner "', fname, '" -c quit() ''',sprintf(' > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt"',exchfolder,simid,exchfolder,simid) ]);
    %         f = fopen(sprintf('%s/sim%d/NeuronLogFile.txt',exchfolder,simid));
    %         txt = fscanf(f,'%c');
    %         fclose(f);
    %         txt
    if ~isempty(oldpwd)
        cd(oldpwd);
    end
end


if ~isempty(strfind(options,'-cl'))
    if ~isempty(answer.StdOut)
        str = regexp(answer.StdOut,'Your job [0-9]*','match','ONCE');
%         ncount = cellfun(@numel,str);
%         [~, ind] = max(ncount);
        jobid = str2double(str(10:end));%{ind}
        % there might be error if many jobs are run, because answer might not
        % be 1
    else
        jobid = NaN;
    end
else
    %     jobid = str2double(regexp(txt,'(?<=ProcessId = )\w*','match'));
    %     if isempty(jobid)
    jobid = NaN;
    %     end
end
end

function minterf = make_nseg(tree,minterf,params,mech)
%does the same as the d_lambda procedure in NEURON
%necessary to find nearest segment which will be calculated
dodlambda = 0;
doeach = 0;

if ischar(params.nseg)
%     pl = Pvec_tree(tree);  % path length of tree..
    len = len_tree(tree);
    idpar = idpar_tree(tree);
    
    if strcmpi(params.nseg,'dlambda')
        dodlambda = 1;
        
        D =  tree.D;
        freq = 100;
        
        if isfield(params,'d_lambda')
            d_lambda = params.d_lambda;
        else
            
            d_lambda = 0.1;
        end
    elseif ~isempty(strfind(params.nseg,'ach')) % check for "each"
        each = double(cell2mat(textscan(params.nseg,'%*s %d'))); % get number
        doeach = 1;
    end
    
end

minterf(:,4) = 0;

for sec = 0:max(minterf(:,2))  %go through all sections
    secstart = find(minterf(:,2) == sec & minterf(:,3) == 0);
    secend = find(minterf(:,2) == sec & minterf(:,3) == 1,1,'last');
    if dodlambda || doeach
        secnodestart = minterf(secstart,1);
        if secnodestart == 0  % this means the current section is the tiny section added for neuron... this should have nseg = 1
            secnodestart = 1;
            flag = true;
        else
            flag = false;
        end
        secnodestart2 = minterf(secstart+1,1);
        secnodeend = minterf(secend,1);
        
        if dodlambda
            if isfield(tree,'rnames') && isfield(mech,tree.rnames{tree.R(secnodestart2)}) && isfield(mech.(tree.rnames{tree.R(secnodestart2)}),'pas') && all(isfield(mech.(tree.rnames{tree.R(secnodestart2)}).pas,{'Ra','cm'}))
                Ra = mech.(tree.rnames{tree.R(secnodestart2)}).pas.Ra;
                cm = mech.(tree.rnames{tree.R(secnodestart2)}).pas.cm;
            elseif isfield(mech,'all') && isfield(mech.all,'pas') && all(isfield(mech.all.pas,{'Ra','cm'}))
                Ra = mech.all.pas.Ra;
                cm = mech.all.pas.cm;
            else
                %NEURON standard values for Ra and cm
                if isfield(tree,'rnames')
                    warning(sprintf('Ra or cm of region %s in tree %s not specified',tree.rnames{tree.R(secnodestart)},tree.name),'Ra or cm not specified','replace')
                else
                    warning('Cannot find passive parameters for nseg calculation! If this is desired, you should define a fixed nseg value','No passive paramers found','replace')
                end
                Ra = 35.4;
                cm = 1;
            end
        end

        if flag
            L = 0.0001;   % this is the length according to root_tree
        else
%             L =  pl(secnodeend) - pl(secnodestart); %length of section
            L = sum(len(secnodestart2:secnodeend)); %length of section
        end
        if dodlambda
%             lambda_f = 0;
%             %from here same calculation as in fixnseg
%             for in = secnodestart2:secnodeend
%                 if in == secnodestart2   % if lastnode was a branching node it is not in a row with next node.account for that
%                     lambda_f = lambda_f + (pl(in)-pl(secnodestart))/sqrt(D(secnodestart)+D(in));
%                 else
%                     lambda_f = lambda_f + (pl(in)-pl(in-1))/sqrt(D(in-1)+D(in));
%                 end
%             end
            lambda_f = sum(len(secnodestart2:secnodeend)./sqrt(D(idpar(secnodestart2:secnodeend)) + D(secnodestart2:secnodeend)));
            
            lambda_f = lambda_f * sqrt(2) * 1e-5*sqrt(4*pi*freq*Ra*cm);
            
            if lambda_f == 0
                lambda_f = 1;
            else
                lambda_f = L/lambda_f;
            end
            %         fprintf('%g\n',lambda_f)
            nseg = floor((L/(d_lambda*lambda_f)+0.9)/2)*2 + 1;     %!%!%! recheck this in NEURON book
        else
            nseg = floor(L/each)+1;
        end
    else
        nseg = params.nseg;
    end
    %     fprintf('%d\n',nseg);
    if isfield(params,'accuracy')
        if params.accuracy == 2 || (params.accuracy == 1 && (~isempty(strfind(tree.rnames(tree.R(minterf(secend))),'axon')) || ~isempty(strfind(tree.rnames(tree.R(minterf(secend))),'soma'))) ) %triple nseg if accuracy is necessary also finds axonh
            nseg = 3 * nseg;
        end
    end
    
    %     fprintf('%d\n',nseg)
    pos = (2 * (1:nseg) - 1) / (2*nseg);    %calculate positions
    fac = (secend-secstart+1)/nseg;
    if fac > 1   % more tree nodes than segments
%         if numel(pos) < 10
%             for in = secstart+1:secend  %!%!%! secstart+1 because first segment gets NaN
%                 [~,ind] = min(abs(minterf(in,3) - pos));   %find position to node which is closest to next segment location
%                 minterf(in,4) = pos(ind);                % hier evt ausnahme für anfang und ende der section (=0)
%             end
%         else
            [~,ind] = min(abs(repmat(minterf(secstart+1:secend,3),1,numel(pos)) - repmat(pos,secend-secstart,1)),[],2);  %find position to node which is closest to next segment location
            minterf(secstart+1:secend,4) = pos(ind);
%         end
    else     % more segments than tree nodes, so one has to add entries in minterf
        dupl = 0;
        for p = 1:numel(pos)
            %             [~,ind] = min(abs(minterf(secstart+1:secend+dupl,3) - pos(p)));   %find position to node which is closest to next segment location %!%!%!
            ind = find(min(unique(abs(minterf(secstart+1:secend+dupl,3) - pos(p)))) == abs(minterf(secstart+1:secend+dupl,3) - pos(p)),1,'last'); % might be more cpu-consuming but attaches next pos in right order
            if minterf(secstart+1+ind-1,4) ~= 0
                minterf = minterf([1:secstart+1+ind-1,secstart+1+ind-1:end],:); %duplicate entry
                minterf(secstart+1+ind,4) = pos(p);  %overwrite duplicated entry with new value
                dupl = dupl +1;
            else
                minterf(secstart+1+ind-1,4) = pos(p);                % hier evt ausnahme für anfang und ende der section (=0)
            end
        end
    end
    minterf(secstart,4) = NaN;
end

end


function n = getref(n,neuron,field)
if isfield(neuron{n},field)
    if strcmp(field,'tree')
        if isnumeric(neuron{n}.tree)        % n is already correct ref
            return
        elseif isempty(strfind(neuron{n}.tree,'sim'))
            n = NaN;
            return
        end
        
        while 1
            if ~isnumeric(neuron{n}.tree)  % ref to a sim
                if str2double(neuron{n}.tree(4:end)) == n   % ref to itsself
                    if n == 1
                        n = [];
                        return
                    else
                        n = NaN;
                        break
                    end
                else                        % ref to another sim
                    n = str2double(neuron{n}.tree(4:end));
                end
            else                        % it found the ref, break
                break
            end
        end
        
    else
        while 1
            if isfield(neuron{n},field)
                if isnumeric(neuron{n}.(field))  % ref to a sim
                    if neuron{n}.(field) == n   % ref to itsself
                        n = NaN;
                        break
                    else                        % ref to another sim
                        n = neuron{n}.(field);
                    end
                else                        % it found the ref, break
                    break
                end
            else  % ref to a sim which not that field
                n = NaN;
                break
            end
        end
    end
else
    n = NaN;
end

end