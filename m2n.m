function [out, origminterf,params,tree] = m2n(tree,params,neuron,options)
% function m2n ("Matlab to Neuron") to generate and execute a .hoc-file in NEURON
% with the parameters in the vector <params> and <neuron>;
% The output-file(s) of the NEURON function are read by cn and transferred
% into the output variable out
% Second and third argument are optional; use {} to leave out the second
% argument.
% 'readyflag' is reserved for checking if neuron has finished
% openNeuron [0]: option to open the command window of NEURON (set to '1';
% useful for debugging)
% don't use 'load_file("nrngui.hoc")' in the NEURON-procedure if you don't
% want the NEURON-menue to
% open (and to get the focus --> extremely anyoing...);
%
% neuron.connect =
% neuron.play =  {node , 'param', {tvec},{vec},continbool};
% neuron.APCount= {node, thresh};
%
% options:
%   -w waitbar
%   -d Debug mode (NEURON is opened and some parameters are set)
%   -q quiet mode -> suppress output
%   -cl cluster mode -> files are prepared to be executed on a cluster.
%   -f  force cluster to run neuron directly without qsub
%                       Automatic run is not executed and output files are only read if copied back
%
% This code is based on an idea of Johannes Kasper, a former group-member
% in the Lab of Hermann Cuntz, Frankfurt.
%
% Copyright by marcel.beining@gmail.com, June 2015

%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


interf_file = 'neuron_runthis.hoc'; % the file which will be written

if ~isfield(params,'path')
    params.path = regexprep(pwd,'\\','/');
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

if ~isempty(strfind(options,'-cl'))
    nrn_path = params.server.clpath;
    if ~isfield(params.server,'walltime') || numel(params.server.walltime) ~=3
        params.server.walltime = [0 30 0];
        warndlg('Server walltime not specified correctly (1 by 3 vector in params.server.walltime). Walltime set to 30 minutes');
    end
else
    nrn_path = params.path;
end
if strcmp(nrn_path(end),'/')
    nrn_path = nrn_path(1:end-1);
end


%% initialize basic variables
noutfiles = 0;
readfiles = cell(0);

%% check other input

if nargin < 2 || isempty(params)
    if debug == 1
        %% individual params structure for debug
        params.openNeuron = true;
        params.nseg = 'dlambda';
        params.tstop = 200;
        params.dt = 0.025;
        params.accuracy = 0;
        params.custom = {};
        params.prerun = false;
        
    else
        params = [];
    end
end
% create some standard parameter set
if nargin < 3 || isempty(neuron)
    if debug ~= 1
        warndlg('No input about what to do were given! Standard test is used')
    end
    neuron.pp{1}.IClamp = struct('node', 1, 'del',100,'dur',50,'amp',5);
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

% check trees are correctly inside one cell
if nargin < 1 || isempty(tree)
    errordlg('No tree specified in input')
    out = m2n_error(out,outoptions);
    return
end
if iscell(tree) && iscell(tree{1})
    tree = tree{1};
elseif isstruct(tree)
    tree = {tree};
end

if isfield(params,'morphfolder')
    morphfolder = fullfile(params.path,params.morphfolder);
else
    errordlg('Please give morphfolder in params.morphfolder');
    out = m2n_error(out,outoptions);
    return
end
morphfolder = regexprep(morphfolder,'\\','/');

if ~isfield(params,'neuronpath')
    params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';  % change neuron path if necessary
end

if strfind(options,'-cl')
    if ~isfield(params,'server')
        errordlg('No access data provided for Cluster server. Please specify in params.server')
        out = m2n_error(out,outoptions);
        return
    else
        if isfield(params.server,'connect')
            
        else
            if exist('ganymed-ssh2-build250/ganymed-ssh2-build250.jar','file')
                javaaddpath(which('ganymed-ssh2-build250/ganymed-ssh2-build250.jar'));
            else
                try
                    sshfrommatlabinstall(1)
                catch
                    errordlg('Could not find the ganymed ssh zip file')
                    out = m2n_error(out,outoptions);
                    return
                end
            end
            params.server.connect = sshfrommatlab(params.server.user,params.server.host,params.server.pw);
        end
        if ~isfield(params.server,'clpath')
            %            params.server.clpath = '~';
            %            warndlg('No Path on the Server specified. Root folder will be used')
            errordlg('No Path on Server specified')
            out = m2n_error(out,outoptions);
            return
        end
    end
end



if isfield(params,'exchfolder')
    exchfolder = fullfile(params.path,params.exchfolder);
    if strfind(options,'-cl')
        nrn_exchfolder = fullfile(params.server.clpath,params.exchfolder);
    else
        nrn_exchfolder = exchfolder;
    end
else
    exchfolder = fullfile(params.path,'m2n_exchange');
    if strfind(options,'-cl')
        nrn_exchfolder = fullfile(params.server.clpath,'m2n_exchange');
    else
        nrn_exchfolder = exchfolder;
    end
    params.exchfolder = 'm2n_exchange';
end
nrn_exchfolder = regexprep(nrn_exchfolder,'\\','/');

if ~isfield(params,'cvode')
    params.cvode = false;
end
if ~isfield(params,'use_local_dt')
    params.use_local_dt = 0;
end
if ~isfield(params,'openNeuron')
    params.openNeuron = false;
end
if ~isfield(params,'nseg') || strcmpi(params.nseg,'d_lambda')
    params.nseg = 'dlambda';
end
if ~isfield(params,'tstart')
    params.tstart = 0;
end
if ~isfield(params,'tstop')
    params.tstop = 200;
end
if ~isfield(params,'dt')
    params.dt = 0.025;
end
if ~isfield(params,'accuracy')
    params.accuracy = 0;
end
if ~isfield(params,'custom')
    params.custom = {};
end
if ~isfield(params,'skiprun')
    params.skiprun = false;
end

% when using Parallel Simulation Run, always rewrite every hoc file!


if exist(params.neuronpath,'file') ~= 2
    errordlg(sprintf('No NEURON software (nrniv.exe) found under "%s"\nPlease give correct path using params.neuronpath',params.neuronpath));
    out = m2n_error(out,outoptions,3);
    return
end

if ~isfield(params,'prerun')
    params.prerun = false;
end

if strfind(options,'-q')
    params.openNeuron = 0;
end

if params.cvode && isnumeric(params.dt)
    warning ('m2n:cvode', 'Dt is set but cvode is active. Dt will be ignored');
end

% create the local and server exchange folder
if exist(exchfolder,'dir') == 0
    mkdir(exchfolder);
end
if strfind(options,'-cl')
    [params.server.connect] = sshfrommatlabissue(params.server.connect,sprintf('rm -rf %s',nrn_exchfolder));
    [params.server.connect] = sshfrommatlabissue(params.server.connect,sprintf('mkdir %s',nrn_exchfolder));
end


% check input

thesetrees = cell(numel(neuron),1);
usestreesof = zeros(numel(neuron),1);
flag = false;
bool = cellfun(@(y) isfield(y,'tree'),neuron);
nempty = cellfun(@isempty,neuron);
if any(nempty)
    errordlg(sprintf('The defined simulation #%d is empty, please check\n',find(nempty)))
    out = m2n_error(out,outoptions);
    return
end
switch sum(bool)
    case 1   % use that treeids defined in that one simulation
        if isnumeric(neuron{bool}.tree)
            thesetrees = repmat({neuron{bool}.tree},numel(neuron),1);
            usestreesof = repmat(find(bool),numel(neuron),1);
        else
            x = getref(1,neuron,'tree');
            if isempty(x) % means it refers to itself (maybe due to usage of m2n_as)..use normal trees..
                thesetrees = repmat({1:numel(tree)},numel(neuron),1);
                usestreesof = ones(numel(neuron),1);
                for n = 1:numel(neuron)
                    neuron{n}.tree = thesetrees{n};
                end
            else
                n = find(bool);
                flag = true;
            end
        end
    case 0      % if no trees are given, use trees that are given to m2n in their order...
        thesetrees = repmat({1:numel(tree)},numel(neuron),1);
        usestreesof = ones(numel(neuron),1);
        for n = 1:numel(neuron)
            neuron{n}.tree = thesetrees{n};
        end
    case numel(neuron)
        for n = 1:numel(neuron)
            x = getref(n,neuron,'tree');
            if ~isnan(x)
                thesetrees{n} = neuron{x}.tree;
                usestreesof(n) = x;
            elseif isempty(x)
                thesetrees = repmat({1:numel(tree)},numel(neuron),1);
                usestreesof = ones(numel(neuron),1);
                for n = 1:numel(neuron)
                    neuron{n}.tree = thesetrees{n};
                end
            else
                flag = true;
                break
            end
        end
    otherwise  % if more than one are given, m2n cannot know which trees you want
        n = find(bool);
        flag = true;
end
if flag
    origminterf = [];
    errordlg(sprintf('Error in neuron{%d}.tree, please check\n',n))
    out = m2n_error(out,outoptions);
    return
end
%
for t = 1:numel(tree)
    if isfield(tree{t},'artificial') && ~isfield(tree{t},'NID')
        tree{t}.NID = strcat('cell_',tree{t}.artificial);           % artificial cells only need one morph hoc file which is named cell_ + the name of the artificial cell..
    end
end
if ~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(morphfolder,strcat(x.NID,'.hoc')),'file'),tree))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet or hoc file is missing! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
%         ind = cellfun(@(x) ~isfield(x,'NID'),tree); % find indices to not transformed trees
        ind = cellfun(@(x) ~isfield(x,'NID'),tree) | ~cellfun(@(x) exist(fullfile(morphfolder,strcat(x.NID,'.hoc')),'file'),tree);
        
%         if ~all(ind)
%             ind2 = find(~ind);
%             ind(ind2(find(cellfun(@(x) exist(fullfile(morphfolder,strcat(x.NID,'.hoc')),'file'),tree(ind2))))) = []; % ignore trees that have the hoc file
%         end
        tree(ind) = m2n_writetrees(params,tree(ind),options);
    else
        out = m2n_error(out,outoptions);
        origminterf = [];
        return
    end
end
origminterf = cell(numel(tree),1);
for t = 1:numel(tree)
    origminterf{t} = load(fullfile(morphfolder,sprintf('%s_minterf.dat',tree{t}.NID)));
end


%% start writing hoc file
for n = 1:numel(neuron)
    if strfind(options,'-d')
        tim = tic;
    end
    for t = 1:numel(tree(thesetrees{n}))
        if ~isfield(tree{thesetrees{n}(t)},'artificial')
            x = getref(n,neuron,'mech');
            minterf{thesetrees{n}(t)} = make_nseg(tree{thesetrees{n}(t)},origminterf{thesetrees{n}(t)},params,neuron{x}.mech{thesetrees{n}(t)});
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
        answer = questdlg(sprintf('Error!\n%s seems to be run by another Matlab instance!\nOverwriting might cause errorneous output!\nIf you are sure that there is no simulation running, we can continue and overwrite. Are you sure? ',fullfile(exchfolder,thisfolder)));
        switch answer
            case 'Yes'
                % iamrunning file is kept and script goes on...
            otherwise
                out.error = 1;
                return
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
        [params.server.connect] = sshfrommatlabissue(params.server.connect,sprintf('mkdir %s/%s',nrn_exchfolder,thisfolder));
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
                fprintf(nfile,sprintf('io = nrn_load_dll("lib_mech/%s")\n',params.nrnmech{c}));
            end
        else
            fprintf(nfile,sprintf('io = nrn_load_dll("lib_mech/%s")\n',params.nrnmech));
        end
    else
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
    if ~isempty(params.custom)
        for c = 1:size(params.custom,1)
            if strcmpi(params.custom{c,2},'start') && exist(fullfile(nrn_path,'lib_custom',params.custom{c,1}),'file')
                fprintf(nfile,sprintf('io = load_file("lib_custom/%s")\n',params.custom{c,1}));
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
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_mech.hoc")\n',nrn_exchfolder,sprintf('sim%d',neuron{n}.mech)) );
        else
            fprintf(nfile,'io = xopen("init_mech.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    
    x = getref(n,neuron,'pp');
    if ~isnan(x)
        fprintf(nfile,'// ***** Place Point Processes *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_pp.hoc")\n',nrn_exchfolder,sprintf('sim%d',neuron{n}.pp)) );
        else
            fprintf(nfile,'io = xopen("init_pp.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    x = getref(n,neuron,'con');
    if ~isnan(x)
        fprintf(nfile,'// ***** Define Connections *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_con.hoc")\n',nrn_exchfolder,sprintf('sim%d',neuron{n}.con)) );
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
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_rec.hoc")\n',nrn_exchfolder,sprintf('sim%d',neuron{n}.record)) );
        else  % else write an own file
            fprintf(nfile,'io = xopen("init_rec.hoc")\n' );
        end
        fprintf(nfile,'\n\n');
    end
    
    
    x = getref(n,neuron,'play');
    if ~isnan(x)
        fprintf(nfile,'// ***** Define vector play sites *****\n');
        if x~=n
            fprintf(nfile,sprintf('io = load_file("%s/%s/init_play.hoc")\n',nrn_exchfolder,sprintf('sim%d',neuron{n}.play)) );
        else
            fprintf(nfile,'io = xopen("init_play.hoc")\n' );
        end
    end
    
    fprintf(nfile,'\n\n');
    fprintf(nfile,'// ***** Last settings *****\n');
    fprintf(nfile,'addsurf_spines()\n');
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
    if ~isempty(params.custom)
        for c = 1:size(params.custom,1)
            if strcmpi(params.custom{c,2},'mid') && exist(fullfile(params.path,'lib_custom',params.custom{c,1}),'file')
                fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,params.custom{c,1}));
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
    if ~isempty(params.custom)
        for c = 1:size(params.custom,1)
            if strcmpi(params.custom{c,2},'end') && exist(fullfile(params.path,'lib_custom',params.custom{c,1}),'file')
                fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,params.custom{c,1}));
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
        fprintf(ofile,'// ***** Load cell morphologiy templates and create artificial cells *****\n');
        for t = thesetrees{n}
            % load templates generated by neuron_template_tree, create one
            % instance of them and add them to the cellList
            
            fprintf(ofile,sprintf('io = xopen("%s//%s.hoc")\n',params.morphfolder,tree{thesetrees{n}(t)}.NID) );
            fprintf(ofile, sprintf('cell = new %s()\n', tree{thesetrees{n}(t)}.NID) );
            if isfield( tree{thesetrees{n}(t)},'params')
                fields = fieldnames( tree{thesetrees{n}(t)}.params);
                for f = 1:numel(fields)
                    fprintf(ofile, sprintf('cell.cell.%s = %g\n',fields{f}, tree{thesetrees{n}(t)}.params.(fields{f})));
                end
            end
            
            fprintf(ofile, 'io = cellList.append(cell)\n');
            
        end
        fprintf(ofile, 'objref cell\n');
        
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define nseg for all cells (called in init_mech) *****\n');
        fprintf(ofile, 'proc make_nseg() {\n');
        fprintf(ofile, 'for CELLINDEX = 0, cellList.count -1 {\n');
        fprintf(ofile, 'if (cellList.o(CELLINDEX).is_artificial == 0) {\n');
        if isfield(params,'nseg') && isnumeric(params.nseg)
            fprintf(ofile, 'forsec cellList.o(CELLINDEX).allreg {\n');
            fprintf(ofile, sprintf('nseg = %f\n}\n}\n}\n}\n',round(params.nseg)) );
            if rem(round(params.nseg),2) == 0
                warndlg('nseg is not odd! Please reconsider nseg');
            end
        elseif isfield(params,'nseg') && strcmpi(params.nseg,'dlambda')
            fprintf(ofile, 'geom_nseg()\n}\n}\n}\n');
        else
            fprintf(ofile, '// No nseg specified!!!\n}\n}\n}\n');
            warndlg('nseg has not been specified (correctly?)! nseg is not set!')
        end
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
            for t = 1:numel(thesetrees{n})
                if numel(neuron{n}.mech) >= t && ~isempty(neuron{n}.mech{t})   && ~isfield(tree{thesetrees{n}(t)},'artificial')    % if a mechanism is defined for this tree
                    if isstruct(neuron{n}.mech{t})          % input must be a structure
                        fields = fieldnames(neuron{n}.mech{t});
                    else
                        continue
                    end
                    
                    if any(strcmpi(fields,'all'))
                        str = sprintf('forsec cellList.o(%d).allreg {\n',t-1);   %neuron:go through this region
                        mechs = fieldnames(neuron{n}.mech{t}.all);                % mechanism names are the fieldnames in the structure
                        for m = 1:numel(mechs)      % loop through mechanisms
                            str = sprintf('%sinsert %s\n',str,mechs{m});        % neuron:insert this mechanism
                            if ~isempty(neuron{n}.mech{t}.all.(mechs{m}))
                                mechpar = fieldnames(neuron{n}.mech{t}.all.(mechs{m}));
                                for p = 1:numel(mechpar)  % loop through mechanism parameters
                                    if strcmpi(mechpar{p},'cm') || strcmpi(mechpar{p},'Ra') || (~isempty(strfind(mechs{m},'_ion')) &&  (numel(mechpar{p}) <= strfind(mechs{m},'_ion') || (numel(mechpar{p}) > strfind(mechs{m},'_ion') && ~strcmp(mechpar{p}(strfind(mechs{m},'_ion')+1),'0'))))       %if mechanism is an ion or passive cm/Ra, leave out mechansim suffix
                                        str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p}));   %neuron: define values
                                    else
                                        str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.all.(mechs{m}).(mechpar{p}));   %neuron: define values
                                    end
                                end
                            end
                        end
                        fprintf(ofile,sprintf('%s}\n\n',str));
                    end
                    
                    if isfield(tree{thesetrees{n}(t)},'R')
                        uR = unique(tree{thesetrees{n}(t)}.R); % Region indices that exist in tree
                        if ~isempty(intersect(tree{thesetrees{n}(t)}.rnames(uR),fields)) %isstruct(neuron{n}.mech{t}.(fields{1}))  %check if mechanism was defined dependent on existent region
                            regs = fields;  %if yes (some of) the input are the regions
                            regs = intersect(tree{thesetrees{n}(t)}.rnames(uR),regs);  % only use those region names which are existent in tree
                            for r = 1 : numel(regs)
                                str = sprintf('forsec cellList.o(%d).reg%s {\n',t-1,regs{r});   %neuron:go through this region
                                mechs = fieldnames(neuron{n}.mech{t}.(regs{r}));                % mechanism names are the fieldnames in the structure
                                for m = 1:numel(mechs)      % loop through mechanisms
                                    str = sprintf('%sinsert %s\n',str,mechs{m});        % neuron:insert this mechanism
                                    
                                    if ~isempty(neuron{n}.mech{t}.(regs{r}).(mechs{m}))
                                        mechpar = fieldnames(neuron{n}.mech{t}.(regs{r}).(mechs{m}));
                                        for p = 1:numel(mechpar)  % loop through mechanism parameters
                                            if strcmpi(mechpar{p},'cm') || strcmpi(mechpar{p},'Ra') || (~isempty(strfind(mechs{m},'_ion')) &&  (numel(mechpar{p}) <= strfind(mechs{m},'_ion') || (numel(mechpar{p}) > strfind(mechs{m},'_ion') && ~strcmp(mechpar{p}(strfind(mechs{m},'_ion')+1),'0'))))       %if mechanism is an ion or passive cm/Ra, leave out mechansim suffix
                                                str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                            else
                                                if numel(neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p})) == 1
                                                    str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                                else
                                                    errordlg(sprintf('Parameter %s of mechanism %s in region %s has more than one value, please check.',mechpar{p},mechs{m},regs{r}))
                                                    out = m2n_error(out,outoptions);
                                                    for nn = 1:numel(neuron)
                                                        delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                                    end
                                                    return
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
                        if ~isfield(tree{thesetrees{n}(t)},'artificial')
                            %                             str = '';
                            [~, ia] = unique(minterf{thesetrees{n}(t)}(:,[2,4]),'rows','stable');  % find real segments in neuron simulation
                            ia = ia(~isnan(minterf{thesetrees{n}(t)}(ia,4))); % remove start nodes of a segment (cause their value belongs to segment -1)
                            ia(numel(ia)+1) = size(minterf{thesetrees{n}(t)},1)+1;   % add one entry
                            
                            mechs = fieldnames(neuron{n}.mech{t}.range);
                            for m = 1:numel(mechs)
                                vars = fieldnames(neuron{n}.mech{t}.range.(mechs{m}));
                                %                                 allvals = zeros(3,0);
                                %                                 thesevars = '';
                                for r = 1:numel(vars)
                                    if numel(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})) == numel(tree{thesetrees{n}(t)}.X)
                                        allvals = zeros(3,0);
                                        for in = 1:numel(ia)-1
                                            thisval = nanmean(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})(minterf{thesetrees{n}(t)}(ia(in),1):minterf{thesetrees{n}(t)}(ia(in+1)-1,1))); % the value is the mean of all tree nodes which are simulated by this segment, if first node is start of section, ignore this one, since it belongs to old region
                                            
                                            if ~isnan(thisval)
                                                allvals = cat(2,allvals,[minterf{thesetrees{n}(t)}(ia(in),[2,4]),thisval]');
                                            end
                                        end
                                        
                                        %                                         thesevars = sprintf('%s"%s_%s",',thesevars,vars{r},mechs{m});
                                        secname = sprintf('range_%s_%s_%s_sec.dat',tree{thesetrees{n}(t)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,secname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(1,:));
                                        fclose(f);
                                        segname = sprintf('range_%s_%s_%s_seg.dat',tree{thesetrees{n}(t)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,segname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(2,:));
                                        fclose(f);
                                        valname = sprintf('range_%s_%s_%s_val.dat',tree{thesetrees{n}(t)}.NID,vars{r},mechs{m});
                                        f = fopen(fullfile(exchfolder,thisfolder,valname) ,'Wt');
                                        fprintf(f,'%g\n',allvals(3,:));
                                        fclose(f);
                                        %                                         fprintf(ofile,sprintf('set_range(%d,"%s","%s","%s","%s_%s")\n',t-1,secname,segname,valname,vars{r},mechs{m}));
                                        rangestr = sprintf('%sset_range(%d,"%s","%s","%s","%s_%s")\n',rangestr,t-1,secname,segname,valname,vars{r},mechs{m});
                                    else
                                        errordlg('Range variable definition should be a vector with same number of elements as tree has nodes')
                                        out = m2n_error(out,outoptions);
                                        for nn = 1:numel(neuron)
                                            delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                                        end
                                        return
                                    end
                                end
                             end
                            
                        else
                            errordlg('Setting range variables for artificial cells is invalid')
                            out = m2n_error(out,outoptions);
                            for nn = 1:numel(neuron)
                                delete(fullfile(exchfolder,sprintf('sim%d',nn),'iamrunning'));   % delete the running mark
                            end
                            return
                        end
                    end
                end
            end
            fprintf(ofile,'\n\n');
            fprintf(ofile,'// ***** Now adjust number of segments *****\n');
            fprintf(ofile,'make_nseg()\n');
            if ~isempty(rangestr)
                fprintf(ofile,'\n\n');
                fprintf(ofile,'// ***** Set specified range variables *****\n');
                fprintf(ofile,rangestr);
                fprintf(ofile,'\n\n');
            end
            if isfield(params,'celsius')
                fprintf(ofile,'\n\nobjref q10\nq10 = new Temperature()\n' ) ;
                fprintf(ofile,sprintf('io = q10.correct(%g)\n\n',params.celsius) ) ;
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
            for t = 1:numel(thesetrees{n})
                if numel(neuron{n}.pp) >= t && ~isempty(neuron{n}.pp{t})   && ~isfield(tree{thesetrees{n}(t)},'artificial')    % if point processes are defined for this tree
                    ppfield = fieldnames(neuron{n}.pp{t});
                    for f1 = 1:numel(ppfield)
                        %%%%
                        node = neuron{n}.pp{t}.(ppfield{f1}).node;
                        for in = 1:numel(node)
                            inode = find(minterf{thesetrees{n}(t)}(:,1) == node(in),1,'first');    %find the index of the node in minterf
                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',t-1,minterf{thesetrees{n}(t)}(inode,2) ) );    % corresponding section of node
                            fprintf(ofile,sprintf('{pp = new %s(%f)\n',ppfield{f1},minterf{thesetrees{n}(t)}(inode,3) ) );  % new pp
                            fields = setdiff(fieldnames(neuron{n}.pp{t}.(ppfield{f1})),{'node','id'});
                            
                            if any(strcmp(ppfield{f1},{'IClamp','SEClamp','SEClamp2','VClamp'})) && (any(strcmp(fields,'times')) || (any(strcmp(fields,'dur')) && numel(neuron{n}.pp{t}.(ppfield{f1}).dur) > 3))
                                if any(strcmp(fields,'times'))
                                    times = sprintf('%f,',neuron{n}.pp{t}.(ppfield{f1}).times);
                                else   %bugfix since seclamp can only use up to 3 duration specifications
                                    times = sprintf('%f,',[0 cumsum(neuron{n}.pp{t}.(ppfield{f1}).dur(1:end-1))]);
                                end
                                amps = sprintf('%f,',neuron{n}.pp{t}.(ppfield{f1}).amp);
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
                                fields = setdiff(fields,{'times','amp','dur'});
                                fprintf(ofile,'io = playtList.append(playt)\n');
                                fprintf(ofile,'io = playList.append(play)\n');
                                fprintf(ofile, 'objref play\n');
                                fprintf(ofile, 'objref playt\n');
                            end
                            
                            for f2 =1:numel(fields)  % go through all parameter fields and declare them
                                if any(strcmpi(fields{f2},{'dur','amp'})) && any(strcmp(ppfield{f1},{'SEClamp2','SEClamp','VClamp'}))   % for dur and amp, there are multiple values
                                    for ff = 1:numel(neuron{n}.pp{t}.(ppfield{f1}).(fields{f2}))
                                        switch ppfield{f1}
                                            case 'VClamp'
                                                fprintf(ofile,sprintf('pp.%s[%d] = %f \n',fields{f2},ff-1,neuron{n}.pp{t}.(ppfield{f1}).(fields{f2})(ff)));
                                            case {'SEClamp','SEClamp2'}
                                                fprintf(ofile,sprintf('pp.%s%d = %f \n',fields{f2},ff,neuron{n}.pp{t}.(ppfield{f1}).(fields{f2})(ff)) );
                                        end
                                    end
                                    
                                else
                                    fprintf(ofile,sprintf('pp.%s = %f \n', fields{f2},neuron{n}.pp{t}.(ppfield{f1}).(fields{f2})) );
                                end
                            end
                            
                            fprintf(ofile,'}\n');
                            fprintf(ofile,'io = ppList.append(pp)\n' );  %append pp to ppList
                            neuron{n}.pp{t}.(ppfield{f1}).id(in) = count;   % save id to pplist in Neuron (for find the correct object for recording later)
                            count = count +1; %ppnum(t) = ppnum(t) +1;  % counter up
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
            for c = 1:size(neuron{n}.con,1)
                str = cell(0);
                nodeflag = false;
                sourcefields = setdiff(fieldnames(neuron{n}.con(c).source),{'cell','watch'});
                
                cell_source = neuron{n}.con(c).source.cell;
                if isempty(sourcefields)   % probably an artificial cell...
                    for t = 1:numel(cell_source)
                        if ~isempty(cell_source(t)) && isfield(tree{thesetrees{n}(cell_source(t))},'artificial')
                            str{t} = sprintf('con = new NetCon(cellList.o(%d).cell,',cell_source(t)-1);
                            
                        else
                            str = sprintf('con = new NetCon(nil,');
                        end
                    end
                else
                    if any(strcmp(sourcefields,'pp'))  % point process is the source
                        x = getref(n,neuron,'pp');
                        pp = neuron{n}.con(c).source.pp;
                        [~,iid] = intersect(neuron{x}.pp{cell_source}.(pp).node,neuron{n}.con(c).source.node); % get reference to the node location of the PPs that should be connected
                        
                        for ii = 1:numel(iid)
                            str{ii} = sprintf('con = new NetCon(ppList.o(%d),',neuron{x}.pp{cell_source}.(pp).id(iid(ii)));
                        end
                    else   % a normal section is the source
                        node = neuron{n}.con(c).source.node;
                        for in = 1:numel(node)
                            inode = find(minterf{thesetrees{n}(cell_source)}(:,1) == node(in),1,'first');    %find the index of the node in minterf
                            if isfield(neuron{n}.con(c).source,'watch') && ischar(neuron{n}.con(c).source.watch)
                                str{in} = sprintf('cellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&%s(%f),',cell_source-1,minterf{thesetrees{n}(cell_source)}(inode,2),neuron{n}.con(c).source.watch,minterf{thesetrees{n}(cell_source)}(inode,3));
                            else
                                str{in} = sprintf('cellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&v(%f),',cell_source-1,minterf{thesetrees{n}(cell_source)}(inode,2),minterf{thesetrees{n}(cell_source)}(inode,3));
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
                                if ~isempty(cell_target(t2)) && isfield(tree{thesetrees{n}(cell_target(t2))},'artificial')
                                    newstr{count} = sprintf('%scellList.o(%d).cell',str{t1},cell_target-1);
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
                        [~,iid] = intersect(neuron{x}.pp{cell_target}.(pp).node,neuron{n}.con(c).target(it).node); % get reference to the node location of the PPs that should be connected
                        for t1 = 1:numel(cell_source)
                            for ii = 1:numel(iid)
                                newstr{count} = sprintf('%sppList.o(%d)',str{t1},neuron{x}.pp{cell_target}.(pp).id(iid(ii)));
                                count = count +1;
                            end
                        end
                        
                    else   % nothing...
                        warndlg('No target specified as connection')
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
            for t = 1:numel(thesetrees{n})
                if numel(neuron{x}.record) >= t && ~isempty(neuron{x}.record{t})  % if a recording site was defined for  this tree
                    recfields = fieldnames(neuron{x}.record{t});
                    if isfield(tree{thesetrees{n}(t)},'artificial')
                        if numel(recfields) > 1 && strcmp(recfields,'record')
                            neuron{x}.record{t} = struct(tree{thesetrees{n}(t)}.artificial,neuron{x}.record{t}); % put the structure in field named as the artificial neuron
                        end
                    end
                    recfields = fieldnames(neuron{x}.record{t});
                    
                    for f1 = 1:numel(recfields)
                        if isfield(tree{thesetrees{n}(t)},'artificial')
                            rectype = 'artificial';
                        elseif strcmp(recfields{f1},'cell') %   size(neuron{x}.record{t},2)<3 || isempty(neuron{x}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                            rectype = 'cell';
                        else
                            rectype = 'pp';
                        end
                        
                        for r = 1:numel(neuron{x}.record{t}.(recfields{f1})) %.record)  % go through all variables to be recorded
                            
                            if strcmp(rectype,'cell')
                                if isfield(tree{thesetrees{n}(t)},'R') && isfield(tree{thesetrees{n}(t)},'rnames')
                                    Rs = tree{thesetrees{n}(t)}.rnames(tree{thesetrees{n}(t)}.R(neuron{x}.record{t}.cell(r).node));       % all region names of trees nodes
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
                                            warndlg(sprintf('Region(s) "%s" of tree %d do not contain mechanism "%s" for recording. Recording in this region is ignored',str(1:end-1),t,strs{end}))
                                        end
                                    end
                                end
                            end
                            
                            if ~any(strcmp(rectype,'artificial'))
                                inode = zeros(numel(neuron{x}.record{t}.(recfields{f1})(r).node),1);
                                for in = 1:numel(neuron{x}.record{t}.(recfields{f1})(r).node)
                                    inode(in) = find(minterf{thesetrees{n}(t)}(:,1) == neuron{x}.record{t}.(recfields{f1})(r).node(in),1,'first');    %find the index of the node in minterf
                                end
                                [realrecs,~,ic] = unique(minterf{thesetrees{n}(t)}(inode,[2,4]),'rows');
                                % put warndlg here !
                            end
                            
                            switch rectype
                                case 'cell'
                                    for in = 1:size(realrecs,1)
                                        %             fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {\n',t-1,minterf{thesetrees{n}(t)}(inode,2) ) );
                                        fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('rec.label("%s at location %06.4f of section %d of cell %d")\n', neuron{x}.record{t}.cell(r).record , realrecs(in,2), realrecs(in,1) ,t-1) ); % label the vector for plotting
                                        if params.cvode
                                            fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&%s(%f),rec,rect)}\n',t-1,realrecs(in,1), neuron{x}.record{t}.cell(r).record, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{x}.record
                                        else
                                            fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),tvec)\n',t-1,realrecs(in,1), neuron{x}.record{t}.cell(r).record, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{x}.record
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
                                            fprintf(ofile,sprintf('rec.label("%s of %s Point Process at location %06.4f of section %d of cell %d")\n', neuron{x}.record{t}.(recfields{f1})(r).record , recfields{f1} , fliplr(realrecs(in,:)) ,t-1) ); % label the vector for plotting
                                            if params.cvode
                                                fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                                fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&ppList.o(%d).%s,rec,rect)}\n',t-1, realrecs(in,1),neuron{x3}.pp{t}.(recfields{f1})(r).id(ind), neuron{x}.record{t}.(recfields{f1})(r).record ) ); % record the parameter x at site y as specified in neuron{x}.record
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
                                            ic(ic == in) = [];  % if node does not correspond to some specified pp at that place, delete it
                                            ic(ic >= in) = ic(ic >= in) - 1;
                                        end
                                        if ~isempty(delin)
                                            realrecs(delin,:) = [];  % if node does not correspond to some specified pp at that place, delete it
                                        end
                                    end
                                    neuron{x}.record{t}.(recfields{f1})(r).rrecs = realrecs;
                                    neuron{x}.record{t}.(recfields{f1})(r).irrecs = ic; % gives the the index to realrecs for each node
                                case 'artificial'
                                    fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                    fprintf(ofile,sprintf('rec.label("%s of artificial cell %s (cell #%d)")\n', neuron{x}.record{t}.cell(r).record , tree{thesetrees{n}(t)}.artificial, t-1) ); % label the vector for plotting
                                    if strcmpi(neuron{x}.record{t}.cell(r).record,'on')
                                        fprintf(ofile,sprintf('nilcon = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',t-1,0.5) );    % for art. cells, make netcon with threshold 0.5
                                        fprintf(ofile,sprintf('io = nilcon.record(rec)\n'));
                                        fprintf(ofile,'io = nilconList.append(nilcon)\n\n' );  %append recording vector to recList
                                        
                                    else
                                        if params.cvode
                                            fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                            fprintf(ofile,sprintf('io = cvode.record(&cellList.o(%d).cell.%s,rec,rect)\n',t-1, neuron{x}.record{t}.cell(r).record ) );  % record the parameter x of artificial cell t-1
                                        else
                                            fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).cell.%s,tvec)\n', t-1, neuron{x}.record{t}.cell.record{r} ) ); % record the parameter x of artificial cell t-1
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
            for t = 1:numel(thesetrees{n})
                if numel(neuron{x2}.APCount) >= t && ~isempty(neuron{x2}.APCount{t})   % if a recording site was defined for  this tree
                    for r = 1: size(neuron{x2}.APCount{t},1)
                        if ~isfield(tree{thesetrees{n}(t)},'artificial')
                            inode = find(minterf{thesetrees{n}(t)}(:,1) == neuron{x2}.APCount{t}(r,1),1,'first');    %find the index of the node in minterf
                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',t-1,minterf{thesetrees{n}(t)}(inode,2) ) );    % corresponding section of node
                            fprintf(ofile,sprintf('{APC = new APCount(%f)\n',minterf{thesetrees{n}(t)}(inode,3) ) );    % make APCCount at position x
                            fprintf(ofile,sprintf('APC.thresh = %f\n',neuron{x2}.APCount{t}(r,2) ) ); % set threshold of APCount [mV]
                        else
                            fprintf(ofile,sprintf('APC = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',t-1,neuron{x2}.APCount{t}(r,2) ) );    % for art. cells, make netcon with threshold
                        end
                        fprintf(ofile,'APCrec = new Vector()\n');
                        fprintf(ofile,'io = APCrecList.append(APCrec)\n');
                        fprintf(ofile,'io = APC.record(APCrecList.o(APCrecList.count()-1))\n');
                        
                        if ~isfield(tree{thesetrees{n}(t)},'artificial')
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
    
    
 %% write init_play.hoc (copied!)
    
    x = getref(n,neuron,'play');
    if (~isnan(x)) && x==n      %rewrite only if one or both of play/APCount are not taken from previous sim or if both are taken from another sim but from different ones (not possible because both are in one hoc)
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_play.hoc') ,'wt');   %open play hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define play sites *****\n');
        %     end
            count = 0;  % counter for playing vector List
            countt = 0; % counter for time vector List
            for t = 1:numel(thesetrees{n})
                if numel(neuron{x}.play) >= t && ~isempty(neuron{x}.play{t})  % if a playing site was defined for  this tree
                    
                    if isfield(tree{thesetrees{n}(t)},'artificial')
                        playfields = fieldnames(neuron{x}.play{t});
                        if numel(playfields) > 1 && strcmp(playfields,'play')
                            neuron{x}.play{t} = struct(tree{thesetrees{n}(t)}.artificial,neuron{x}.play{t}); % put the structure in field named as the artificial neuron
                        end
                    end
                    playfields = fieldnames(neuron{x}.play{t});
                    
                    for f1 = 1:numel(playfields)
                        if isfield(tree{thesetrees{n}(t)},'artificial')
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
                                if isfield(tree{thesetrees{n}(t)},'R') && isfield(tree{thesetrees{n}(t)},'rnames')
                                    Rs = tree{thesetrees{n}(t)}.rnames(tree{thesetrees{n}(t)}.R(neuron{x}.play{t}.cell(r).node));       % all region names of trees nodes
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
                                            warndlg(sprintf('Region(s) "%s" of tree %d do not contain mechanism "%s" for playing. Playing in this region is ignored',str(1:end-1),t,strs{end}))
                                        end
                                    end
                                end
                            end
                            
                            if ~any(strcmp(playtype,'artificial'))
                                inode = zeros(numel(neuron{x}.play{t}.(playfields{f1})(r).node),1);
                                for in = 1:numel(neuron{x}.play{t}.(playfields{f1})(r).node)
                                    inode(in) = find(minterf{thesetrees{n}(t)}(:,1) == neuron{x}.play{t}.(playfields{f1})(r).node(in),1,'first');    %find the index of the node in minterf
                                end
                                [realplays,~,ic] = unique(minterf{thesetrees{n}(t)}(inode,[2,4]),'rows');
                                % put warndlg here !
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
                                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_at_%d_cell_%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,t-1)),'w');
                                        fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).times(1:end-1));
                                        fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).times(end));
                                        fclose(f);
                                        fprintf(ofile,'f = new File()\n');
                                        fprintf(ofile,sprintf('f.ropen("plt_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in),t-1));  %vector file is opened
                                        fprintf(ofile,'playt.scanf(f)\n');    % file is read into time vector
                                        fprintf(ofile,'io = f.close()\n');     %file is closed
                                        fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append playing time vector to playtList
                                        
                                        fprintf(ofile,sprintf('play = new Vector(%f)\n',length(neuron{n}.play{t}.(playfields{f1})(r).value) ) );    % create new playing vector
                                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('pl_%s_at_%d_cell_%d.dat', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,t-1)),'w');
                                        fprintf(f,'%g ', neuron{n}.play{t}.(playfields{f1})(r).value(1:end-1));
                                        fprintf(f,'%g\n', neuron{n}.play{t}.(playfields{f1})(r).value(end));
                                        fclose(f);
                                        fprintf(ofile,'f = new File()\n');
                                        fprintf(ofile,sprintf('f.ropen("pl_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}.(playfields{f1})(r).play , ic(in) ,t-1));  %vector file is opened
                                        fprintf(ofile,'play.scanf(f)\n');     % file is read into play vector
                                        fprintf(ofile,'io = f.close()\n');   %file is closed
                                        fprintf(ofile,sprintf('play.label("playing %s at node %d of cell %d")\n', neuron{n}.play{t}.(playfields{f1})(r).play  , ic(in) ,t-1) ); % label the vector for plotting
                                        fprintf(ofile,sprintf('play.play(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),playtList.o(playtList.count()-1),%d)\n',t-1,realplays(in,1), neuron{n}.play{t}.(playfields{f1})(r).play, realplays(in,2), cont ) ); % play the parameter x at site y as specified in neuron{n}.play
                                        fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList

                                    end
                                    neuron{x}.play{t}.cell(r).rplays = realplays; % gives the section and segment to the playings
                                    neuron{x}.play{t}.cell(r).irplays = ic; % gives the the index to realplays for each node
                                case 'pp'
                                    delin = [];
                                    for in =  1:size(realplays,1)%numel(neuron{x}.play{t}{r,1})  % CAUTION might be wrong
                                        ind = find(neuron{x3}.pp{t}.(playfields{f1}).node == neuron{x}.play{t}.(playfields{f1})(r).node(in));
                                        if ~isempty(ind)
                                            fprintf(ofile,sprintf('play = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new playing vector
                                            fprintf(ofile,sprintf('play.label("%s of %s Point Process at location %06.4f of section %d of cell %d")\n', neuron{x}.play{t}.(playfields{f1})(r).play , playfields{f1} , fliplr(realplays(in,:)) ,t-1) ); % label the vector for plotting
                                            if params.cvode
                                                fprintf(ofile,sprintf('playt = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new playing vector
                                                fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.play(&ppList.o(%d).%s,play,playt)}\n',t-1, realplays(in,1),neuron{x3}.pp{t}.(playfields{f1})(r).id(ind), neuron{x}.play{t}.(playfields{f1})(r).play ) ); % play the parameter x at site y as specified in neuron{x}.play
                                            else
                                                fprintf(ofile,sprintf('io = play.play(&ppList.o(%d).%s,tvec)\n',neuron{x3}.pp{t}.(playfields{f1})(r).id(ind), neuron{x}.play{t}.(playfields{f1})(r).play ) ); % play the parameter x at site y as specified in neuron{x}.play
                                            end
                                            fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
                                            if params.cvode
                                                fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append time playing vector to playList
                                                neuron{x}.play{t}.(playfields{f1})(r).idt(in) = countt;   % reference to find playing in playList
                                                countt = countt +1;
                                            end
                                            neuron{x}.play{t}.(playfields{f1})(r).id(in) = count;   % reference to find playing in playList
                                            count = count +1;
                                        else
                                            delin = cat(1,delin,in);
                                            ic(ic == in) = [];  % if node does not correspond to some specified pp at that place, delete it
                                            ic(ic >= in) = ic(ic >= in) - 1;
                                        end
                                        if ~isempty(delin)
                                            realplays(delin,:) = [];  % if node does not correspond to some specified pp at that place, delete it
                                        end
                                    end
                                    neuron{x}.play{t}.(playfields{f1})(r).rplays = realplays;
                                    neuron{x}.play{t}.(playfields{f1})(r).irplays = ic; % gives the the index to realplays for each node
                                case 'artificial'
                                    fprintf(ofile,sprintf('play = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new playing vector
                                    fprintf(ofile,sprintf('play.label("%s of artificial cell %s (cell #%d)")\n', neuron{x}.play{t}.cell(r).play , tree{thesetrees{n}(t)}.artificial, t-1) ); % label the vector for plotting
                                    if strcmpi(neuron{x}.play{t}.cell(r).play,'on')
                                        fprintf(ofile,sprintf('nilcon = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',t-1,0.5) );    % for art. cells, make netcon with threshold 0.5
                                        fprintf(ofile,sprintf('io = nilcon.play(play)\n'));
                                        fprintf(ofile,'io = nilconList.append(nilcon)\n\n' );  %append playing vector to playList
                                        
                                    else
                                        if params.cvode
                                            fprintf(ofile,sprintf('playt = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new playing vector
                                            fprintf(ofile,sprintf('io = cvode.play(&cellList.o(%d).cell.%s,play,playt)\n',t-1, neuron{x}.play{t}.cell(r).play ) );  % play the parameter x of artificial cell t-1
                                        else
                                            fprintf(ofile,sprintf('io = play.play(&cellList.o(%d).cell.%s,tvec)\n', t-1, neuron{x}.play{t}.cell.play{r} ) ); % play the parameter x of artificial cell t-1
                                        end
                                        if params.cvode
                                            fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append time playing vector to playList
                                            neuron{x}.play{t}.cell(r).idt = countt;   % reference to find playing in playList
                                            countt = countt +1;
                                        end
                                    end
                                    fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
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
    %     if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
    ofile = fopen(fullfile(exchfolder,thisfolder,'save_rec.hoc') ,'wt');   %open record hoc file in write modus
    fprintf(ofile,'// * Write Recordings to Files *\n');
    %     end
    x = getref(n,neuron,'record');
    x2 = getref(n,neuron,'APCount');
    x3 = getref(n,neuron,'pp');
    if ~isnan(x)
        out{n}.record = cell(1,numel(thesetrees{n}));   % initialize output of cn
        
        for t = 1:numel(thesetrees{n})
            if numel(neuron{x}.record) >= t && ~isempty(neuron{x}.record{t})
                %                 for r = 1: size(neuron{x}.record{t},1)
                recfields = fieldnames(neuron{x}.record{t});
                
                for f1 = 1:numel(recfields)
                    if isfield(tree{thesetrees{n}(t)},'artificial')
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
                                    fname = sprintf('cell%d_sec%d_loc%06.4f_%s', t-1, neuron{x}.record{t}.cell(r).rrecs(in,:), neuron{x}.record{t}.(recfields{f1})(r).record  );
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
                                    %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
                                    %                                     readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{x}.record{t}.(recfields{f1}).record{r} , neuron{x}.record{t}.(recfields{f1}).rrecs(in,1), neuron{x}.record{t}.(recfields{f1}).rrecs(in,2) };
                                    readfiles{noutfiles} = {sprintf('%s.dat',fname), n, t , 'cell', neuron{x}.record{t}.(recfields{f1})(r).record , neuron{x}.record{t}.(recfields{f1})(r).node(neuron{x}.record{t}.(recfields{f1})(r).irrecs == in) }; %neuron{x}.record{t}.(recfields{f1})(r).node(in) };
                                end
                            case 'pp'
                                for in = 1:size(neuron{x}.record{t}.(recfields{f1})(r).rrecs,1)
                                    fname = sprintf('%s_cell%d_sec%d_loc%06.4f_%s',recfields{f1},t-1, neuron{x}.record{t}.(recfields{f1})(r).rrecs(in,:), neuron{x}.record{t}.(recfields{f1})(r).record  );
                                    
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
                                    %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
                                    %                                     readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{x}.record{t}.(recfields{f1}).record{r} };
                                    readfiles{noutfiles} = {sprintf('%s.dat',fname), n, t , recfields{f1}, neuron{x}.record{t}.(recfields{f1})(r).record , neuron{x}.record{t}.(recfields{f1})(r).node(neuron{x}.record{t}.(recfields{f1})(r).irrecs == in)};%neuron{x}.record{t}.(recfields{f1})(r).node(in) };
                                end
                            case 'artificial'
                                fname = sprintf('cell%d_%s',t-1, neuron{x}.record{t}.(recfields{f1})(r).record );
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
                                %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{x}.record{t}{r,1},neuron{x}.record{t}{r,2} ) , 'record' ,  t , neuron{x}.record{t}{r,2} ,neuron{x}.record{t}{r,1} };
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
        for t = 1:numel(thesetrees{n})
            if numel(neuron{x}.APCount) >= t && ~isempty(neuron{x}.APCount{t})     % if a recording site was defined for  this tree
                for r = 1: size(neuron{x}.APCount{t},1)
                    fname = sprintf('cell%d_node%d_APCtimes.dat',t-1,neuron{x}.APCount{t}(r,1) );
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
    
    
    if strfind(options,'-cl') %transfer files to server
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
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{7});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{7});
            m = m + 1;
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{8});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{8});
            m = m + 1;
        end
        if getref(n,neuron,'play') == n
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{9});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{9});
            m = m + 1;
        end
        if isempty(strfind(options,'-f'))
            %create job
            ofile = fopen(fullfile(exchfolder,thisfolder,'start_nrn.job') ,'wt');
            
            fprintf(ofile,'#!/bin/bash\n');
            fprintf(ofile,'# write standard output to file\n');
            fprintf(ofile,sprintf('#PBS -o simstart_%s.oe\n',regexprep(datestr(now),{' ','\-','\:'},'_')));
            fprintf(ofile,'# calculate for 30 minutes on 5 core, max. 512 MB of RAM per process\n');
            fprintf(ofile,sprintf('#PBS -l walltime=%s,nodes=5,pmem=512m\n',sprintf('%02d:%02d:%02d',params.server.walltime)));
            fprintf(ofile,'# load needed modules \n');
            fprintf(ofile,'module load openmpi/gcc/64/1.3.3\n');
            fprintf(ofile,'# change to path with your executale\n');
            %     fprintf(ofile,sprintf('cd %s\n',nrn_exchfolder));
            fprintf(ofile,'# start your program with mpirun with 5 processes\n');
            fprintf(ofile,sprintf('mpirun -np 5 nrngui -nobanner -nogui -mpi %s//%s \n',nrn_exchfolder,interf_file));
            fclose(ofile);
            localfilename{m} = fullfile(exchfolder,'start_nrn.job');
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,'start_nrn.job');
        end
        params.server.connect = sftpfrommatlab(params.server.connect,localfilename,remotefilename);
    end
    
    if strfind(options,'-d')
        tim = toc(tim);
        fprintf(sprintf('Sim %d: HOC writing time: %g min %.2f sec\n',n,floor(tim/60),rem(tim,60)))
    end
end

%% Execute NEURON
if ~isempty(strfind(options,'-cl'))
    warndlg('Multi-Thread and server execution not yet implemented')
    out = m2n_error(out,outoptions);
    return
else
    num_cores = feature('numCores');
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
            if isempty(strfind(options,'-f'))
                pause(0.8);
                for s = find(simids==1)
                    pause(0.2);
                    [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('qstat %d',jobid(s)));
                    
                    if ~isempty(answer{1})      %job not finished yet
                        answer = textscan(answer{3},'%*[^QR]%s%*[^QR]');
                        if strcmpi(answer{1},'R')
                            if ~flag(s)
                                display(sprintf('Simulation %d is calculated on cluster',s))
                                if strfind(options,'-d')
                                    tim = toc(tim);
                                    fprintf(sprintf('Cluster queue wait time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
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
                % with direct mpirun matlab waits for command to be finished so there is no need
                % for while loop..
                oanswer = answer;
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
                        errordlg(sprintf('There was an error in Simulation %d:\n******************************\n%s\n******************************\nDue to that m2n has no output to that Simulation.',s(ss),txt(1:min(numel(txt),2000))));
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
                answer = questdlg(sprintf('Waitbar was closed, m2n stopped continuing. Only finished data is returned. If accidently, retry.\nClose all NEURON instances?\n (Caution if several Matlab instances are running)'),'Close NEURON instances?','Close','Ignore','Ignore');
                if strcmp(answer,'Close')
                   system('taskkill /F /IM nrniv.exe');
                end
%                 errordlg('Waitbar was closed, m2n stopped continuing. Only finished data is returned. If accidently, retry.')
                simids(simids<2) = 4;
%                 out = m2n_error(out,outoptions);
                fclose all;
                
            end
        end
    end
    
    if strfind(options,'-cl')
        s = find(simids==1);
        timm = tic;
        for ss = 1:numel(s)
            [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/%s/readyflag',nrn_exchfolder,sprintf('sim%d',s(ss))));
            if isempty(answer)    % then there was an error during executing NEURON
                if strfind(options,'-f')
                    errordlg(sprintf('There was an error during NEURON simulation:\n %s.',strcat(oanswer{:})))
                else
                    [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('ls *e%d*',jobid(s(ss))));
                    if ~isempty(answer)
                        % was planned to directly show error but...
                        %                 scptomatlab(params.server.connect,exchfolder,answer{1})
                        %                 f = fopen(fullfile(exchfolder,answer{1}));
                        %                 errfile = textscan(f,'%s','Delimiter','\n');
                        %                 errstr =
                        errordlg(sprintf('There was an error during NEURON simulation. Please refer to cluster output file "%s".',answer{1}))
                    end
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
        elseif simids(readfiles{f}{2}) == 4  % m2n was aborted 
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
%                 errordlg('Waitbar was closed during data loading. If accidently, retry.')
                out = m2n_error(out,outoptions);
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
    if (outoptions.nocell && ~isfield(out,'error')) || ~outoptions.nocell && any(cellfun(@(x) ~isfield(x,'error'),out))
        out = m2n_error(out,outoptions,0);
    end
    if strfind(options,'-d')
        tim = toc(tim);
        fprintf(sprintf('Data loading time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
end

for n = 1:numel(neuron)
    delete(fullfile(exchfolder,sprintf('sim%d',n),'iamrunning'));   % delete the running mark
end


end


function [jobid,tim] = exec_neuron(simid,exchfolder,nrn_exchfolder,interf_file,params,options)
%% Execute NEURON
if strfind(options,'-d')
    tim = tic;
else
    tim=[];
end
switch params.openNeuron
    case 1
        opsign = '&';%sprintf(' > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt" &',exchfolder,simid,exchfolder,simid); % execute the program in foreground and hand control to NEURON
    case 0
        opsign = sprintf(' -c quit() > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt" \n exit',exchfolder,simid,exchfolder,simid);  % execute the program iconified
end
% execute the file in neuron:
fname = regexprep(fullfile(exchfolder,sprintf('sim%d',simid),interf_file),'\\','/');

if ~isempty(strfind(options,'-q'))
    if ~isempty(strfind(options,'-cl'))
        if strfind(options,'-f')
            [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('nrngui -nobanner -nogui %s/%s/%s \n',nrn_exchfolder,sprintf('sim%d',simid),interf_file));  %!%!%!%!
        else
            [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('qsub %s/%s/%s',nrn_exchfolder,sprintf('sim%d',simid),'start_nrn.job'));
            fprintf(sprintf('Answer server after submitting: %s\nExtracing Job Id and wait..\n',answer{1}))
        end
    else
        oldpwd = '';
        if ~strcmpi(pwd,params.path)
            oldpwd = pwd;
            cd(params.path);   % in order to have NEURON the path as its starting folder
        end
        system(['start /B ' params.neuronpath ' -nobanner "' fname sprintf('" -c quit() > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt"',exchfolder,simid,exchfolder,simid)]); %&,char(13),'exit&']); %nrniv statt neuron
%         system(['wmic process call create ''', params.neuronpath, ' -nobanner "', fname, '" -c quit() ''',sprintf(' > "%s/sim%d/NeuronLogFile.txt" 2> "%s/sim%d/ErrorLogFile.txt"',exchfolder,simid,exchfolder,simid) ]);
%         f = fopen(sprintf('%s/sim%d/NeuronLogFile.txt',exchfolder,simid));
%         txt = fscanf(f,'%c');
%         fclose(f); 
%         txt
        if ~isempty(oldpwd)
            cd(oldpwd);
        end
    end
else
    if (strfind(options,'-cl'))
        % this is also quiet since no neuron gui was installed...
        if strfind(options,'-f')
            [params.server.connect,answer] = sshfrommatlabissue(params.server.connect,sprintf('mpirun -np 5 nrngui -nobanner -nogui -mpi %s/%s \n',nrn_exchfolder,interf_file));
            
        else
            [params.server.connect,answer] = shfrommatlabissue(params.server.connect,sprintf('qsub %s/%s',nrn_exchfolder,'start_nrn.job'));
            fprintf(sprintf('Answer server after submitting: %s\n',answer{1}))
        end
    else
        oldpwd = '';
        if ~strcmpi(pwd,params.path)
            oldpwd = pwd;
            cd(params.path);
        end
        dos([params.neuronpath ' -nobanner "' fname '" ' opsign]); %&,char(13),'exit&']); %nrniv statt neuron
        if ~isempty(oldpwd)
            cd(oldpwd);
        end
    end
end
if ~isempty(strfind(options,'-cl'))
    if numel(answer) == 1 && isempty(strfind(options,'-f'))
        str = regexp(answer{1},'[0-9]*','match');
        ncount = cellfun(@numel,str);
        [~, ind] = max(ncount);
        jobid = str2double(str{ind});
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
if ischar(params.nseg) && strcmpi(params.nseg,'dlambda')
    dodlambda = 1;
    pl = Pvec_tree(tree);  % path length of tree..
    D =  tree.D;
    freq = 100;
    
    if isfield(params,'d_lambda')
        d_lambda = params.d_lambda;
    else
        
        d_lambda = 0.1;
    end
else
    dodlambda = 0;
end

minterf(:,4) = 0;

for sec = 0:max(minterf(:,2))  %go through all sections
    
    secstart = find(minterf(:,2) == sec & minterf(:,3) == 0);
    secend = find(minterf(:,2) == sec & minterf(:,3) == 1,1,'last');
    
    if dodlambda
        secnodestart = minterf(secstart,1);
        if secnodestart == 0  % this means the current section is the tiny section added for neuron... this should have nseg = 1
            secnodestart = 1;
            flag = true;
        else
            flag = false;
        end
        secnodestart2 = minterf(secstart+1,1);
        if isfield(tree,'rnames') && isfield(mech,tree.rnames{tree.R(secnodestart2)}) && isfield(mech.(tree.rnames{tree.R(secnodestart2)}),'pas') && all(isfield(mech.(tree.rnames{tree.R(secnodestart2)}).pas,{'Ra','cm'}))
            Ra = mech.(tree.rnames{tree.R(secnodestart2)}).pas.Ra;
            cm = mech.(tree.rnames{tree.R(secnodestart2)}).pas.cm;
        elseif isfield(mech,'all') && isfield(mech.all,'pas') && all(isfield(mech.all.pas,{'Ra','cm'}))
            Ra = mech.all.pas.Ra;
            cm = mech.all.pas.cm;
        else
            %NEURON standard values for Ra and cm
            if isfield(tree,'rnames')
                warndlg(sprintf('Ra or cm of region %s in tree %s not specified',tree.rnames{tree.R(secnodestart)},tree.name),'Ra or cm not specified','replace')
            else
                warndlg('Cannot find passive parameters for nseg calculation! If this is desired, you should define a fixed nseg value','No passive paramers found','replace')
            end
            Ra = 35.4;
            cm = 1;
        end
        secnodeend = minterf(secend,1);
        if flag
            L = 0.0001;   % this is the length according to root_tree
        else
            L =  pl(secnodeend) - pl(secnodestart); %length of section
        end
        lambda_f = 0;
        %from here same calculation as in fixnseg
        for in = secnodestart2:secnodeend
            if in == secnodestart2   % if lastnode was a branching node it is not in a row with next node.account for that
                lambda_f = lambda_f + (pl(in)-pl(secnodestart))/sqrt(D(secnodestart)+D(in));
            else
                lambda_f = lambda_f + (pl(in)-pl(in-1))/sqrt(D(in-1)+D(in));
            end
        end
        lambda_f = lambda_f * sqrt(2) * 1e-5*sqrt(4*pi*freq*Ra*cm);
        
        if lambda_f == 0
            lambda_f = 1;
        else
            lambda_f = L/lambda_f;
        end
        %         fprintf('%g\n',lambda_f)
        nseg = floor((L/(d_lambda*lambda_f)+0.9)/2)*2 + 1;     %!%!%! recheck this in NEURON book
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
        for in = secstart+1:secend  %!%!%! secstart+1 because first segment gets NaN
            [~,ind] = min(abs(minterf(in,3) - pos));   %find position to node which is closest to next segment location
            minterf(in,4) = pos(ind);                % hier evt ausnahme für anfang und ende der section (=0)
        end
        
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
        end
    end
else
    n = NaN;
end

end
