function [out, minterf,params,tree] = m2n(tree,params,neuron,options)
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
% in the morphological modeling lab Frankfurt.
%
% Copyright by marcel.beining@gmail.com, April 2015

% check trees are correctly inside one cell
if nargin < 1 || isempty(tree)
    errordlg('No tree specified in input')
    return
end
if iscell(tree) && iscell(tree{1})
    tree = tree{1};
elseif isstruct(tree)
    tree = {tree};
end

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
numcell = numel(tree);
noutfiles = 0;
readfiles = cell(0);
orderchanged = false;
changed = struct('morph',0,'stim',1,'basic',1,'lib',1,'rec',1,'play',1,'mech',1,'pp',1,'con',1);


%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


interf_file = 'neuron_runthis.hoc'; % the file which will be written

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
        params.changed = changed;
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
    neuron.stim{1} = {'IClamp', 1, stuct('del',100,'dur',50,'amp',5)};
    neuron.record{1} = {1 , 'v' , 'node'};
    neuron.APCount{1} = {1,-30}; % {node, tresh}
    neuron.mech.all.pas = [];
    neuron.mech.soma.hh = [];
    neuron.mech.axon.hh = [];
end

if isstruct(neuron)         % transform structure neuron to cell neuron
    if numel(neuron) == 1
        neuron = {neuron};
    else
        neuron = arrayfun(@(x) {x},neuron);
    end
end
out = cell(numel(neuron),1);

if ~isfield(params,'neuronpath')
    params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';  % change neuron path if necessary
end

if strfind(options,'-cl')
    if ~isfield(params,'server')
        errordlg('No access data provided for Cluster server. Please specify in params.server')
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
                    return
                end
            end
            params.server.connect = sshfrommatlab(params.server.user,params.server.host,params.server.pw);
        end
        if ~isfield(params.server,'clpath')
            %            params.server.clpath = '~';
            %            warndlg('No Path on the Server specified. Root folder will be used')
            errordlg('No Path on Server specified')
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
if ~isfield(params,'changed')
    params.changed = changed;
elseif ~isempty(setdiff(fieldnames(changed),fieldnames(params.changed)))
    fld = fieldnames(changed);
    for f = 1:numel(fld)
        if ~isfield(params.changed,fld{f})
            params.changed.(fld{f}) = 1;
        end
    end
end
% when using Parallel Simulation Run, always rewrite every hoc file!
if numel(neuron) > 1
    fld = fieldnames(params.changed);
    for f = 1:numel(fld)
        params.changed.(fld{f}) = 1;
    end
end

if exist(params.neuronpath,'file') ~= 2
    errordlg(sprintf('No NEURON software (nrniv.exe) found under "%s"\nPlease give correct path using params.neuronpath',params.neuronpath));
    for n = 1:numel(neuron)
        out{n}.error = 3;
    end
    return
end


if isfield(params,'morphfolder')
    morphfolder = fullfile(params.path,params.morphfolder);
else
    morphfolder = exchfolder;
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
    params.changed.morph = 1;
elseif ~exist(fullfile(exchfolder,'init_cells.hoc'),'file')
    params.changed.morph = 1;
end
if ~isfield(params,'prerun')
    params.prerun = false;
end
if ~isfield(params,'access')
    params.access = [find(~cellfun(@(x) isfield(x,'artificial'),tree),1,'first') 1];      % std accessing first non-artificial tree at node 1
end


if strfind(options,'-q')
    params.openNeuron = 0;
elseif strfind(options,'-d')
    params.changed.basic = 1;
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

badchars = 0;
minterf = cell(numcell,1);
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
        tname = params.tname;
        tflag = true;
    elseif artflag && ~isfield(tree{t},'name')
        tname = tree{t}.artificial;
        tflag = true;
    elseif isfield(tree{t},'name')
        tname = tree{t}.name;
        tflag = false;
    else
        tname = 'Tree';
        tflag = true;
    end
    if any(strfind(tname,'%'))
        badchars = badchars +numel(strfind(tname,'%'));
        tname(strfind(tname,'%')) = [];
    end
    if any(strfind(tname,'.'))
        badchars = badchars +numel(strfind(tname,'.'));
        tname(strfind(tname,'.')) = '_';
    end
    
    tname = strcat('cell_',tname);
    if tflag
        tname = sprintf('%s_%d%d',tname,floor(t/10),rem(t,10));
    end
    for n = 1:numel(neuron)
        neuron{n}.cellIDs{t} = tname;  % save tree names as cellIDs for neuron
    end
    if strfind(options,'-cl')
        [params.server.connect, answer] = sshfrommatlabissue(params.server.connect,sprintf('ls %s/%s.hoc',nrn_morphfolder,tname));
        fchk =  ~isempty(answer{1});
    else
        fchk = exist(fullfile(morphfolder,sprintf('%s.hoc',tname)),'file');
    end
    
    if params.changed.morph || fchk == 0     % if morphology does not already exists
        oname = tname;
        [tname, nix, minterf{t}] = neuron_template_tree (tree{t}, fullfile(morphfolder,sprintf('%s.hoc',tname)), [], '-m');
        if strfind(options,'-cl')   %transfer files to server
            params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s.hoc',oname)),sprintf('%s/%s.hoc',nrn_morphfolder,oname));
            pause(0.1)
            params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.dat',oname)),sprintf('%s/%s_minterf.dat',nrn_morphfolder,oname));
            pause(0.1)
            params.server.connect = sftpfrommatlab(params.server.connect,fullfile(morphfolder,sprintf('%s_minterf.mat',oname)),sprintf('%s/%s_minterf.mat',nrn_morphfolder,oname));
        end
    else
        minterf{t} = load(fullfile(morphfolder,sprintf('%s_minterf.dat',tname)));
    end
    if ~artflag
        Ra = 0;
        cm = 0;
        for n = 1:numel(neuron)     % just try to find the parameter set with maximum Ra and cm levels and adjust nseg to that
            try
                Ra = max(Ra,neuron{n}.mech{t}.all.pas.Ra);
                cm = max(cm,neuron{n}.mech{t}.all.pas.cm);
            end
            fnames = fieldnames(neuron{n}.mech{t});
            for f = 1:numel(fnames)
                try
                    Ra = max(Ra,neuron{n}.mech{t}.(fnames{f}).pas.Ra);
                    cm = max(cm,neuron{n}.mech{t}.(fnames{f}).pas.cm);
                end                
            end
        end
        if cm > 0 && Ra > 0
            mech.all.pas.Ra = Ra;
            mech.all.pas.cm = cm;
            minterf{t} = make_nseg(tree{t},minterf{t},params,mech);
        else
            minterf{t} = make_nseg(tree{t},minterf{t},params,[]);
        end
    end
end
if badchars > 0
    %     warndlg(sprintf('Caution! %d bad chars had to be removed or replaced from the tree names since they cause writing errors! Please be sure to not use "%%" and "." in the names',badchars),'Bad characters removed');
end
if strfind(options,'-d')
    tim = toc(tim);
    fprintf(sprintf('Tree writing time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
end

%% define some output parameters...commented out, instead load tvec from neuron to have exact ts (also necessary for cvode)
% % out.t = params.tstart:params.dt:params.tstop;
% % this is unfortunately necessary to let Matlab have the same tvec size as
% % NEURON, since NEURON produces roundoff errors during long runs
% out.t = NaN(params.tstop/params.dt+1,1);
% out.t(1) = params.tstart;
% ind = 2;
% flag=false;
% while 1
%     if out.t(ind-1) + params.dt <= params.tstop
%         out.t(ind) = out.t(ind-1) + params.dt;
%         ind = ind +1;
%     else
%         if isnan(out.t(end))    % this is due roundoff error -> one entry less..
%             out.t(end) = [];
%             flag = true;
%         end
%         break
%     end
% end
% if isempty(strfind(options,'-r'))       % make nice numbers if not intended to use correct t values
%     out.t = params.tstart:params.dt:params.tstop-(flag*params.dt);
% end

%% start writing hoc file
for n = 1:numel(neuron)
    thisfolder = sprintf('sim%d',n);
    
    % delete the readyflag if it exists
    if exist(fullfile(exchfolder,thisfolder),'dir') == 0
        mkdir(fullfile(exchfolder,thisfolder));
    end
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
    if strfind(options,'-d')
        tim = tic;
    end
    if params.changed.basic || params.changed.lib || params.changed.morph     %rewrite only if something has changed influencing this file
        nfile = fopen(fullfile(exchfolder,thisfolder,interf_file) ,'wt');   %open resulting hoc file in write modus
        
        fprintf(nfile,'// ***** This is a NEURON hoc file automatically created by the Matlab-NEURON interface. *****\n');
        fprintf(nfile,'// ***** Copyright by Marcel Beining and Johannes Kasper, Clinical Neuroanatomy, Goethe University Frankfurt*****\n\n');
        %initialize variables in NEURON
        fprintf(nfile,'// General variables: i, CELLINDEX, debug_mode, accuracy\n');
        
        fprintf(nfile,'// ***** Initialize Variables *****\n');
        fprintf(nfile,'objref f\n');
        fprintf(nfile,'objref nil,cvode,strf,tvec,cell,cellList,pp,ppList,stim,stimList,con,conList,nilcon,nilconList,rec,recList,rect,rectList,playt,playtList,play,playList,APCrec,APCrecList,APC,APCList,APCcon,APCconList \n cellList = new List() // comprises all instances of cell templates, also artificial cells\n stimList = new List() // comprises all SEClamp/VClamp/IClamp electrodes introduced into any cell\n ppList = new List() // comprises all Point Processes of any cell\n conList = new List() // comprises all NetCon objects\n recList = new List() //comprises all recording vectors\n rectList = new List() //comprises all time vectors of recordings\n playtList = new List() //comprises all time vectors for play objects\n playList = new List() //comprises all vectors played into an object\n APCList = new List() //comprises all APC objects\n APCrecList = new List()\n APCconList = new List() //comprises all APC recording vectors\n cvode = new CVode() //the Cvode object\n');%[',numel(tree),']\n'  ;
        
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
            fprintf(nfile,sprintf('io = f.wopen("%s//%s//tvec.dat")\n',nrn_exchfolder,thisfolder)  );  % open file for this time vector with write perm.
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
        fprintf(nfile, sprintf('io = xopen("%s/lib_genroutines/fixnseg.hoc")\n',nrn_path) );
        fprintf(nfile, sprintf('io = xopen("%s/lib_genroutines/genroutines.hoc")\n',nrn_path) );
        fprintf(nfile, sprintf('io = xopen("%s/lib_genroutines/pasroutines.hoc")\n',nrn_path) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Load custom libraries *****\n');
        if ~isempty(params.custom)
            for c = 1:size(params.custom,1)
                if strcmpi(params.custom{c,2},'start') && exist(fullfile(nrn_path,'lib_custom',params.custom{c,1}),'file')
                    fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,params.custom{c,1}));
                end
            end
        end
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Load cell morphologies and create artificial cells *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_cells.hoc")\n',nrn_exchfolder,thisfolder) );
        %     fprintf(nfile,'\n\n');
        %     fprintf(nfile,'// ***** Load passive model *****\n');
        %     fprintf(nfile,sprintf('io = xopen("%s/init_pas.hoc")\n',nrn_exchfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Load mechanisms *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_mech.hoc")\n',nrn_exchfolder,thisfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Place Point Processes *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_pp.hoc")\n',nrn_exchfolder,thisfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Define Connections *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_con.hoc")\n',nrn_exchfolder,thisfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Place stimulations *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_stim.hoc")\n',nrn_exchfolder,thisfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Define recording sites *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_rec.hoc")\n',nrn_exchfolder,thisfolder) );
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Define vector play sites *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/init_play.hoc")\n',nrn_exchfolder,thisfolder) );
        
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
        if numel(params.access) > 1 % if there is any non-artificial cell defined
            fprintf(nfile,sprintf('access cellList.o(%d).allregobj.o(%d).sec\n',params.access(1)-1,minterf{1}(params.access(2),2)) );
        end
        fprintf(nfile,'make_nseg()\n');
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// ***** Include prerun or standard run replacing custom code *****\n');
        if ~isempty(params.custom)
            for c = 1:size(params.custom,1)
                if strcmpi(params.custom{c,2},'mid') && exist(fullfile(params.path,'lib_custom',params.custom{c,1}),'file')
                    fprintf(nfile,sprintf('io = load_file("%s/lib_custom/%s")\n',nrn_path,params.custom{c,1}));
                    %                 if size(params.custom(c),2) > 2 && strcmpi(params.custom{c,3},'skiprun')
                    %                     skiprun = true;
                    %                 end
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
        fprintf(nfile,'// ***** Write Data to Files *****\n');
        fprintf(nfile,sprintf('io = xopen("%s/%s/save_rec.hoc")\n',nrn_exchfolder,thisfolder) );
        
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
        fprintf(nfile,sprintf('io = f.wopen("%s/%s/%s")\n',nrn_exchfolder,thisfolder,'readyflag' ) );       % create the readyflag file
        fprintf(nfile,'io = f.close()\n');   % close the filehandle
        if ~params.openNeuron
            fprintf(nfile,'quit()\n');  % exit NEURON if it was defined so in the parameters
        end
        
        fprintf(nfile,'\n\n');
        fprintf(nfile,'// *-*-*-*-* END *-*-*-*-*\n');
        
        fclose(nfile);
    end
    %% write init_cells.hoc
    if params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_cells.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Load cell morphologiy templates and create artificial cells *****\n');
        for t = 1:numel(tree)
            % load templates generated by neuron_template_tree, create one
            % instance of them and add them to the cellList
            
            fprintf(ofile,sprintf('io = xopen("%s//%s.hoc")\n',nrn_morphfolder,neuron{n}.cellIDs{t}) );
            fprintf(ofile, sprintf('cell = new %s()\n', neuron{n}.cellIDs{t}) );
            fprintf(ofile, 'io = cellList.append(cell)\n');
            
        end
        fprintf(ofile, 'objref cell\n');
        
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define nseg for all cells *****\n');
        fprintf(ofile, 'proc make_nseg() {\n');
        fprintf(ofile, 'for CELLINDEX = 0, cellList.count -1 {\n');
        fprintf(ofile, 'if (cellList.o(CELLINDEX).is_artificial == 0) {\n');
        %     fprintf(ofile, 'if (strf.is_artificial(cellList.o(CELLINDEX)) == 0) {\n');
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
    end

    if params.changed.mech || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_mech.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Insert mechanisms *****\n');
        if isfield(neuron{n},'mech')
            for t = 1:numel(tree)
                if numel(neuron{n}.mech) >= t && ~isempty(neuron{n}.mech{t})   && ~isfield(tree{t},'artificial')    % if a mechanism is defined for this tree
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
                            %                         if size(neuron{n}.mech{t}.all.(mechs{m}),2) == 2     %if parameter definition is no mby2 cell array,rearrange
                            %                             mechpars = neuron{n}.mech{t}.(regs{r}).(mechs{m});
                            %                         else
                            %                             resh = size(neuron{n}.mech{t}.all.(mechs{m}),2)/2;
                            %                             mechpars = reshape(neuron{n}.mech{t}.all.(mechs{m}),2,resh)';
                            %                         end
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
                    
                    if isfield(tree{t},'R')
                        uR = unique(tree{t}.R); % Region indices that exist in tree
                        if ~isempty(intersect(tree{t}.rnames(uR),fields)) %isstruct(neuron{n}.mech{t}.(fields{1}))  %check if mechanism was defined dependent on existent region
                            regs = fields;  %if yes (some of) the input are the regions
                            regs = intersect(tree{t}.rnames(uR),regs);  % only use those region names which are existent in tree
                            for r = 1 : numel(regs)
                                str = sprintf('forsec cellList.o(%d).reg%s {\n',t-1,regs{r});   %neuron:go through this region
                                mechs = fieldnames(neuron{n}.mech{t}.(regs{r}));                % mechanism names are the fieldnames in the structure
                                for m = 1:numel(mechs)      % loop through mechanisms
                                    str = sprintf('%sinsert %s\n',str,mechs{m});        % neuron:insert this mechanism
                                    
                                    %                             if size(neuron{n}.mech{t}.(regs{r}).(mechs{m}),2) == 2     %if parameter definition is no mby2 cell array,rearrange
                                    %                                 mechpars = neuron{n}.mech{t}.(regs{r}).(mechs{m});
                                    %                             else
                                    %                                 resh = size(neuron{n}.mech{t}.(regs{r}).(mechs{m}),2)/2;
                                    %                                 mechpars = reshape(neuron{n}.mech{t}.(regs{r}).(mechs{m}),2,resh)';
                                    %                             end
                                    if ~isempty(neuron{n}.mech{t}.(regs{r}).(mechs{m}))
                                        mechpar = fieldnames(neuron{n}.mech{t}.(regs{r}).(mechs{m}));
                                        for p = 1:numel(mechpar)  % loop through mechanism parameters
                                            if strcmpi(mechpar{p},'cm') || strcmpi(mechpar{p},'Ra') || (~isempty(strfind(mechs{m},'_ion')) &&  (numel(mechpar{p}) <= strfind(mechs{m},'_ion') || (numel(mechpar{p}) > strfind(mechs{m},'_ion') && ~strcmp(mechpar{p}(strfind(mechs{m},'_ion')+1),'0'))))       %if mechanism is an ion or passive cm/Ra, leave out mechansim suffix
                                                str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                            else
                                                str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.(regs{r}).(mechs{m}).(mechpar{p}));   %neuron: define values
                                            end
                                        end
                                    end
                                end
                                fprintf(ofile,sprintf('%s}\n\n',str));
                            end
                        end
                    end
                    
                    if any(strcmpi(fields,'range')) 
                        if ~isfield(tree{t},'artificial')
                            fprintf(ofile,'\n\n');
                            fprintf(ofile,'// ***** Set specified range variables *****\n');
                            str = '';
                            [nix, ia] = unique(minterf{t}(:,[2,4]),'rows','stable');
                            ia(end+1) = numel(tree{t}.X)+1;
                            
                            mechs = fieldnames(neuron{n}.mech{t}.range);
                            
                            for in = 1:numel(ia)-1
                                if in == 1 || minterf{t}(ia(in),2) ~= minterf{t}(ia(in-1),2)
                                    substr = sprintf('cellList.o(%d).allregobj.o(%d).sec {\n',t-1,minterf{t}(ia(in),2));
                                    strflag = false;
                                end
                                
                                for m = 1:numel(mechs)
                                    vars = fieldnames(neuron{n}.mech{t}.range.(mechs{m}));
                                    for r = 1:numel(vars)
                                        if any(numel(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})) == [numel(tree{t}.X) 1])
                                            if numel(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})) == 1
                                                flag = true;
                                            else
                                                flag = false;
                                            end
                                            
                                            if minterf{t}(ia(in),1) == 0 || flag   % accounting for added zero length section and if only one value given
                                                thisval = neuron{n}.mech{t}.range.(mechs{m}).(vars{r})(1);
                                            else
                                                thisval = nanmean(neuron{n}.mech{t}.range.(mechs{m}).(vars{r})(minterf{t}(ia(in),1):minterf{t}(ia(in+1)-1,1)));
                                            end
                                            if ~isnan(thisval)
                                                strflag = true;
                                                substr = sprintf('%s%s_%s(%f) = %f\n',substr,vars{r},mechs{m},minterf{t}(ia(in),3),thisval);
                                            end
                                            
                                            
                                        else
                                            errordlg('Range variable definition should be a vector with same number of elements as tree has nodes')
                                            return
                                        end
                                    end
                                end
                                if minterf{t}(ia(in),2) ~= minterf{t}(ia(in+1),2)
                                    if strflag
                                        str = sprintf('%s%s}\n\n',str,substr);
                                    end
                                    %  str = strcat(str,'\n');
                                end
                            end
                        else
                            errordlg('Setting range variables for artificial cells is invalid')
                            return
                        end
                        fprintf(ofile,sprintf('%s\n\n',str));
                        
                    end
                    %                 % now go through the rest
                    %                 mechs = fields;
                    %                 mechs = setdiff(mechs,tree{t}.rnames); %delete the Region definitions
                    %                 mechs = setdiff(mechs,'all'); %delete the Region definitions   % everything remaining should be unused regions..
                    %                 %CAUTION if regions have been defined which are not in the
                    %                 %tree
                    %                 if ~isempty(mechs)
                    %                     flag=false;
                    %                     for m = 1:numel(mechs) %loop through all (putative) mechanisms
                    %                         if isstruct(neuron{n}.mech{t}.(mechs{m}))
                    %                             continue %workaround to avoid not existent regions, everything else must be a mechanism!
                    %                         end
                    %                         if ~flag   %if it's the first time the loop works,write the forsec argument in neuron to go through all sections
                    %                             str = sprintf('forsec cellList.o(%d).allreg {\n',t-1);
                    %                             flag = true;
                    %                         end
                    %                         str = sprintf('%sinsert %s\n',str,mechs{m});   %neuron: insert this mechanism
                    %
                    % %                         if size(neuron{n}.mech{t}.(mechs{m}),2) == 2     %if parameter definition is no mby2 cell array,rearrange
                    % %                             mechpars = neuron{n}.mech{t}.(mechs{m});
                    % %                         else
                    % %                             resh = size(neuron{n}.mech{t}.(mechs{m}),2)/2;
                    % %                             mechpars = reshape(neuron{n}.mech{t}.(mechs{m}),resh,2)';
                    % %                         end
                    %                         if ~isempty(neuron{n}.mech{t}.(mechs{m}))
                    %                             mechpar = fieldnames(neuron{n}.mech{t}.(mechs{m}));
                    %                             for p = 1:numel(mechpar)         %loop through all mechnism parameters
                    %                                 if strfind(mechs{m},'_ion')        %if mechanism is an ion, leave out mechanism suffix
                    %                                     str = sprintf('%s%s = %g\n',str,mechpar{p},neuron{n}.mech{t}.(mechs{m}).(mechpar{p}));
                    %                                 else
                    %                                     str = sprintf('%s%s_%s = %g\n',str,mechpar{p},mechs{m},neuron{n}.mech{t}.(mechs{m}).(mechpar{p}));
                    %                                 end
                    %                             end
                    %                         end
                    %                     end
                    %                     fprintf(ofile,sprintf('%s}\n\n',str));      %write the string and close section loop
                    %                 end
                end
            end
            if isfield(params,'celsius')
                fprintf(ofile,'\n\nobjref q10\nq10 = new Temperature()\n' ) ;
                fprintf(ofile,sprintf('io = q10.correct(%g)\n\n',params.celsius) ) ;
            end
        end
        fclose(ofile);          %close file
    end
    
    if params.changed.pp || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_pp.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Place synapses or other point processes *****\n');
        if isfield(neuron{n},'pp')
            ppnum = zeros(numel(tree),1);
            for t = 1:numel(tree)
                if numel(neuron{n}.pp) >= t && ~isempty(neuron{n}.pp{t})   && ~isfield(tree{t},'artificial')    % if point processes are defined for this tree
                    % this adds PPs specified...maybe a function
                    % distributing them will do it better later?
                    for s = 1:size(neuron{n}.pp{t},1)
                        inode = find(minterf{t}(:,1) == neuron{n}.pp{t}{s,2},1,'first');    %find the index of the node in minterf
                        fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',t-1,minterf{t}(inode,2) ) );    % corresponding section of node
                        fprintf(ofile,sprintf('{pp = new %s(%f)\n',neuron{n}.pp{t}{s,1},minterf{t}(inode,3) ) );  % new pp
                        if numel(neuron{n}.pp{t}) > 2 && isstruct(neuron{n}.pp{t}{s,3})          % input must be a structure
                            fields = fieldnames(neuron{n}.pp{t}{s,3});
                            for f =1:numel(fields)
                                fprintf(ofile,sprintf('pp.%s = %g\n',fields{f},neuron{n}.pp{t}{s,3}.(fields{f})));
                            end
                        end
                        fprintf(ofile,'}\n');
                        fprintf(ofile,'io = ppList.append(pp)\n' );  %append pp to ppList
                        ppnum(t) = ppnum(t) +1;
                    end
                end
            end
            fprintf(ofile, 'objref pp\n');
        end
        fclose(ofile);
    else            % ppnum has to be redefined anyways
        if isfield(neuron{n},'pp')
            ppnum = zeros(numel(tree),1);
            for t = 1:numel(tree)
                if numel(neuron{n}.pp) >= t && ~isempty(neuron{n}.pp{t})   && ~isfield(tree{t},'artificial')    % if PPs are defined for this tree
                    for s = 1:size(neuron{n}.pp{t},1)
                        ppnum(t) = ppnum(t) +1;
                    end
                end
            end
        end
    end
    
    
    if params.changed.con || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_con.hoc') ,'wt');   %open morph hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define Connections *****\n');
        if isfield(neuron{n},'con')
            % should have look like: {source(node or point process), what to
            % watch, target, threshold, delay, weight}
            for c = 1:size(neuron{n}.con,1)
                str = '';
                nodeflag = false;
                switch neuron{n}.con{c,1}
                    
                    case 'cell'
                        if ischar(neuron{n}.con{c,2})
                            t = str2double(neuron{n}.con{c,2});
                        else
                            t = neuron{n}.con{c,2};
                        end
                        if ~isempty(t) && isfield(tree{t},'artificial')
                            %                         if ~isempty(neuron{n}.con{c,3}) && ischar(neuron{n}.con{c,3})
                            %                             str = sprintf('%scon = new NetCon(&cellList.o(%d)(%s),',str,t-1,neuron{n}.con{c,3});
                            %                         else
                            str = sprintf('%scon = new NetCon(cellList.o(%d).cell,',str,t-1);
                            %                         end
                        else
                            str = sprintf('%scon = new NetCon(nil,',str);
                        end
                    case 'node'
                        expr = regexp(neuron{n}.con{c,2},'\.','split');
                        t = str2double(expr{1});
                        inode = str2double(expr{2});
                        inode = find(minterf{t}(:,1) == inode,1,'first');    %find the index of the node in minterf
                        if ~isempty(neuron{n}.con{c,3}) && ~isnan(neuron{n}.con{c,3})
                            str = sprintf('%scellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&%s(%f),',str,t-1,minterf{t}(inode,2),neuron{n}.con{c,3},minterf{t}(inode,3));
                        else
                            str = sprintf('%scellList.o(%d).allregobj.o(%d).sec {con = new NetCon(&v(%f),',str,t-1,minterf{t}(inode,2),minterf{t}(inode,3));
                        end
                        nodeflag = true;
                    case 'pp'
                        expr = regexp(neuron{n}.con{c,2},'\.','split');
                        t = str2double(expr{1});
                        s = str2double(expr{2});
                        %                     if ~isempty(neuron{n}.con{c,3}) && ischar(neuron{n}.con{c,3})
                        %                         str = sprintf('%scon = new NetCon(&ppList.o(%d)(%s),',str,sum(ppnum(1:t-1))+s-1,neuron{n}.con{c,3});
                        %                     else
                        str = sprintf('%scon = new NetCon(ppList.o(%d),',str,sum(ppnum(1:t-1))+s-1);
                        %                     end
                        
                    otherwise
                        str = sprintf('%scon = new NetCon(nil,',str);
                end
                
                switch neuron{n}.con{c,4}
                    case 'cell'
                        if ischar(neuron{n}.con{c,2})
                            t = str2double(neuron{n}.con{c,5});
                        else
                            t = neuron{n}.con{c,5};
                        end
                        if ~isempty(t) && isfield(tree{t},'artificial')
                            str = sprintf('%scellList.o(%d).cell',str,t-1);
                        else
                            str = sprintf('%snil',str);
                        end
                    case 'pp'
                        expr = regexp(neuron{n}.con{c,5},'\.','split');
                        t = str2double(expr{1});
                        s = str2double(expr{2});
                        str = sprintf('%sppList.o(%d)',str,sum(ppnum(1:t-1))+s-1);
                    otherwise
                        str = sprintf('%snil',str);
                        
                end
                if size(neuron{n}.con(c,:),2) >= 8 && numel(cat(1,neuron{n}.con{c,6:8})) == 3
                    str = sprintf('%s,%g,%g,%g)\n',str,neuron{n}.con{c,6},neuron{n}.con{c,7},neuron{n}.con{c,8});   %threshold , delay,weight
                else
                    str = sprintf('%s)\n',str);
                end
                
                str = sprintf('%sio = conList.append(con)',str);  %append con to conList
                if nodeflag
                    str = sprintf('%s}\n',str);
                else
                    str = sprintf('%s\n',str);
                end
                fprintf(ofile,str);  % new connection
                
                
                
            end
            fprintf(ofile, 'objref con\n');
            %         for t = 1:numel(tree)
            %             if numel(neuron{n}.con) >= t && ~isempty(neuron{n}.con{t})     % if a connection is defined for this tree
            %             end
            %         end
        end
        fprintf(ofile,'\n\n');
        fclose(ofile);
    end
    
    
    %% write init_stim.hoc
    if params.changed.stim || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_stim.hoc') ,'wt');   %open stim hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Place stimulations *****\n');
        if isfield(neuron{n},'stim')
            stimnum = zeros(numel(tree),1);
            for t = 1:numel(tree)
                if numel(neuron{n}.stim) >= t && ~isempty(neuron{n}.stim{t}) && ~isfield(tree{t},'artificial')   % if a stimulation is defined for this tree
                    for s = 1: size(neuron{n}.stim{t},1)
                        inode = find(minterf{t}(:,1) == neuron{n}.stim{t}{s,2},1,'first');    %find the index of the node in minterf
                        fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',t-1,minterf{t}(inode,2) ) );    % corresponding section of node
                        switch neuron{n}.stim{t}{s,1}      % if SEClamp and VClamp dur and amp would be handled equally this could be simplified much more =/
                            case 'IClamp'
                                fprintf(ofile,sprintf('{stim = new IClamp(%f)\n',minterf{t}(inode,4) ) );  % new stim
                                fields = fieldnames(neuron{n}.stim{t}{s,3});
                                if any(strcmp(fields,'times'))
                                    times = sprintf('%f,',neuron{n}.stim{t}{s,3}.times);
                                    times = times(1:end-1);
                                    amps = sprintf('%f,',neuron{n}.stim{t}{s,3}.amp);
                                    amps = amps(1:end-1);
                                    
                                    fprintf(ofile,'playt = new Vector()\n');
                                    fprintf(ofile,sprintf('playt = playt.append(%s)\n',times));
                                    fprintf(ofile,'play = new Vector()\n');
                                    fprintf(ofile,sprintf('play = play.append(%s)\n',amps));
                                    fprintf(ofile,'play.play(&stim.amp,playt)\n');
                                    fprintf(ofile,'stim.dur = 1e15\n');
                                    fprintf(ofile,'stim.del = -1e4\n');
                                    
                                    fprintf(ofile,'io = playtList.append(playt)\n');
                                    fprintf(ofile,'io = playList.append(play)\n');
                                    fprintf(ofile, 'objref play\n');
                                    fprintf(ofile, 'objref playt\n');
                                    %                                 fields = setdiff(fields,{'times','amp','dur'});
                                else
                                    for f = 1:numel(fields)
                                        fprintf(ofile,sprintf('stim.%s = %f \n',fields{f},neuron{n}.stim{t}{s,3}.(fields{f})));
                                    end
                                end
                                %                             fprintf(ofile,sprintf('stim.del = %f \nstim.dur = %f \nstim.amp = %f }\n',neuron{n}.stim{t}{s,3:4} ) );    %define stim properties as described in neuron{n}.stim
                            case 'VClamp'
                                fprintf(ofile,sprintf('{stim = new VClamp(%f)\n',minterf{t}(inode,4) ) );  % new stim
                                fields = fieldnames(neuron{n}.stim{t}{s,3});   % get parameter names
                                if any(strcmp(fields,'times'))
                                    times = sprintf('%f,',neuron{n}.stim{t}{s,3}.times);
                                    times = times(1:end-1);     %delete last comma
                                    amps = sprintf('%f,',neuron{n}.stim{t}{s,3}.amp);
                                    amps = amps(1:end-1);    %delete last comma
                                    
                                    fprintf(ofile,'playt = new Vector()\n');
                                    fprintf(ofile,sprintf('playt = playt.append(%s)\n',times));
                                    fprintf(ofile,'play = new Vector()\n');
                                    fprintf(ofile,sprintf('play = play.append(%s)\n',amps));
                                    fprintf(ofile,'play.play(&stim.amp[0],playt)\n');
                                    fprintf(ofile,'stim.dur[0] = 1e15\n');
                                    
                                    fprintf(ofile,'io = playtList.append(playt)\n');
                                    fprintf(ofile,'io = playList.append(play)\n');
                                    fprintf(ofile, 'objref play\n');
                                    fprintf(ofile, 'objref playt\n');
                                    fields = setdiff(fields,{'times','amp','dur'});
                                end
                                
                                for f = 1:numel(fields)     % loop through all parameters and write them in hoc
                                    if any(strcmpi(fields{f},{'dur','amp'}))    % for dur and amp, there are multiple values
                                        for ff = 1:numel(neuron{n}.stim{t}{s,3}.(fields{f}))
                                            fprintf(ofile,sprintf('stim.%s[%d] = %f \n',fields{f},ff-1,neuron{n}.stim{t}{s,3}.(fields{f})(ff)));
                                        end
                                    else
                                        fprintf(ofile,sprintf('stim.%s = %f \n',fields{f},neuron{n}.stim{t}{s,3}.(fields{f})));
                                    end
                                end
                                
                            case 'SEClamp'
                                fprintf(ofile,sprintf('{stim = new SEClamp(%f)\n',minterf{t}(inode,4) ) );  % new stim
                                fields = fieldnames(neuron{n}.stim{t}{s,3});   % get parameter names
                                if any(strcmp(fields,'times')) || (any(strcmp(fields,'dur')) && numel(neuron{n}.stim{t}{s,3}.dur) > 3)
                                    if any(strcmp(fields,'times'))
                                        times = sprintf('%f,',neuron{n}.stim{t}{s,3}.times);
                                    else   %bugfix since seclamp can only use up to 3 duration specifications
                                        times = sprintf('%f,',[0 cumsum(neuron{n}.stim{t}{s,3}.dur(1:end-1))]);
                                    end
                                    times = times(1:end-1);     %delete last comma
                                    amps = sprintf('%f,',neuron{n}.stim{t}{s,3}.amp);
                                    amps = amps(1:end-1);    %delete last comma
                                    
                                    fprintf(ofile,'playt = new Vector()\n');
                                    fprintf(ofile,sprintf('playt = playt.append(%s)\n',times));
                                    fprintf(ofile,'play = new Vector()\n');
                                    fprintf(ofile,sprintf('play = play.append(%s)\n',amps));
                                    fprintf(ofile,'play.play(&stim.amp1,playt)\n');
                                    fprintf(ofile,'stim.dur1 = 1e15\n');
                                    
                                    fprintf(ofile,'io = playtList.append(playt)\n');
                                    fprintf(ofile,'io = playList.append(play)\n');
                                    fprintf(ofile, 'objref play\n');
                                    fprintf(ofile, 'objref playt\n');
                                    fields = setdiff(fields,{'times','amp','dur'});
                                end
                                
                                for f = 1:numel(fields)     % loop through all parameters and write them in hoc
                                    if any(strcmpi(fields{f},{'dur','amp'}))    % for dur and amp, there are multiple values
                                        for ff = 1:numel(neuron{n}.stim{t}{s,3}.(fields{f}))
                                            fprintf(ofile,sprintf('stim.%s%d = %f \n',fields{f},ff,neuron{n}.stim{t}{s,3}.(fields{f})(ff)));
                                        end
                                    else
                                        fprintf(ofile,sprintf('stim.%s = %f \n',fields{f},neuron{n}.stim{t}{s,3}.(fields{f})));
                                    end
                                end
                            otherwise
                                errordlg(sprintf('Stimulation electrode type %s not know. Please use IClamp,VClamp or SEClamp or define your electrode via neuron.pp',neuron{n}.stim{t}{s,1}))
                                
                        end
                        fprintf(ofile,'} \n io = stimList.append(stim)\n' );  %append stim to stimList
                        stimnum(t) = stimnum(t) +1;
                    end
                    fprintf(ofile,'\n');
                end
                
            end
            fprintf(ofile, 'objref stim\n');
        end
        fclose(ofile);
    else            % stimnum has to be redefined anyways
        if isfield(neuron{n},'stim')
            stimnum = zeros(numel(tree),1);
            for t = 1:numel(tree)
                if numel(neuron{n}.stim) >= t && ~isempty(neuron{n}.stim{t})   && ~isfield(tree{t},'artificial')    % if PPs are defined for this tree
                    for s = 1:size(neuron{n}.stim{t},1)
                        stimnum(t) = stimnum(t) +1;
                    end
                end
            end
        end
    end
    %% write init_rec.hoc
    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_rec.hoc') ,'wt');   %open record hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define recording sites *****\n');
    end
    if isfield(neuron{n},'record')
        for t = 1:numel(tree)
            if numel(neuron{n}.record) >= t && ~isempty(neuron{n}.record{t})  % if a recording site was defined for  this tree
                for r = 1: size(neuron{n}.record{t},1)
                    neuron{n}.record{t}{r,1} = reshape(neuron{n}.record{t}{r,1},numel(neuron{n}.record{t}{r,1}),1);     % put all vectors along first dimension...
                end
                [C,nix,ic] = unique(neuron{n}.record{t}(:,2));
                if numel(C) < numel(neuron{n}.record{t}(:,2))  % multiple assignment of same variable
                    newlin = cell(0,3);
                    delc = false(numel(ic),1);
                    for c = 1:numel(C)
                       if sum(ic==c) > 1                        % there is a multiple assignment
                           newlin = cat(1,newlin,{unique(cat(1,neuron{n}.record{t}{ic==c,1})), C{c}  , neuron{n}.record{t}{find(ic==c,1,'first'),3}  });  % store a new line where all nodes of this recording were put in one
                           delc(ic==c) = true;
                       end                        
                    end  
                    neuron{n}.record{t}(delc,:) = [];           % delete all multiple assignments
                    neuron{n}.record{t} = cat(1,neuron{n}.record{t},newlin);   % add the all-in-one lines
                end
                for r = 1: size(neuron{n}.record{t},1)
                    if isfield(tree{t},'artificial')
                        rectype = 'artificial';
                        if size(neuron{n}.record{t},2) < 2 || isempty(neuron{n}.record{t}{r,2})
                            neuron{n}.record{t}{r,2} = neuron{n}.record{t}{r,1};                    % if parameter to record was defined in the first entry
                        end
                    elseif size(neuron{n}.record{t},2)<3 || isempty(neuron{n}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                        rectype = 'node';
                    else
                        rectype = neuron{n}.record{t}{r,3};
                        if any(strcmp(rectype,{'PP','Pp','syn','Syn'}))
                            rectype = 'pp';
                        end
                        if any(strcmp(rectype,{'Stim','STIM'}))
                            rectype = 'stim';
                        end
                    end
                    if strcmp(rectype,'node')
                        if isfield(tree{t},'R') && isfield(tree{t},'rnames')
                            Rs = tree{t}.rnames(tree{t}.R(neuron{n}.record{t}{r,1}));       % all region names of trees nodes
                            strs = regexp(neuron{n}.record{t}{r,2},'_','split');            % split record string to get mechanism name
                            if numel(strs)>1   %any(strcmp(strs{1},{'v','i'}))             % check if record variable is variable of a mechanism or maybe global
                                ignorethese = false(1,numel(neuron{n}.record{t}{r,1}));
                                uRs = unique(Rs);
                                str = '';
                                for u = 1:numel(uRs)                                            % go through regions to be recorded
                                    if (~isfield(neuron{n}.mech{t},uRs{u}) || isfield(neuron{n}.mech{t},uRs{u}) && ~isfield(neuron{n}.mech{t}.(uRs{u}),strs{end})) &&  (~isfield(neuron{n}.mech{t},'all') || isfield(neuron{n}.mech{t},'all') && ~isfield(neuron{n}.mech{t}.all,strs{end}))             % check if this region also has the mechanism to be recorded
                                        ignorethese = ignorethese | strcmp(uRs{u},Rs);           % if not ignore these region for recording
                                        str = strcat(str,uRs{u},'/');
                                    end
                                end
                                if ~isempty(str)
                                    neuron{n}.record{t}{r,1}(ignorethese) = [];                     % delete the recording nodes which should be ignored
                                    warndlg(sprintf('Region(s) "%s" of tree %d do not contain mechanism "%s" for recording. Recording in this region is ignored',str(1:end-1),t,strs{end}))
                                end
                            end
                        end
                    end
                    
                    
                    inode = zeros(numel(neuron{n}.record{t}{r,1}),1);
                    for in = 1:numel(neuron{n}.record{t}{r,1})
                        inode(in) = find(minterf{t}(:,1) == neuron{n}.record{t}{r,1}(in),1,'first');    %find the index of the node in minterf
                    end
                    realrecs = unique(minterf{t}(inode,[2,4]),'rows');
                    
                    switch rectype
                        case 'node'
                            if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                for in = 1:size(realrecs,1)
                                    %             fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {\n',t-1,minterf{t}(inode,2) ) );
                                    fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                    fprintf(ofile,sprintf('rec.label("%s at location %06.4f of section %d of cell %d")\n', neuron{n}.record{t}{r,2} , realrecs(in,2), realrecs(in,1) ,t-1) ); % label the vector for plotting
                                    if params.cvode
                                        fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&%s(%f),rec,rect)}\n',t-1,realrecs(in,1), neuron{n}.record{t}{r,2}, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{n}.record
                                    else
                                        fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),tvec)\n',t-1,realrecs(in,1), neuron{n}.record{t}{r,2}, realrecs(in,2) ) ); % record the parameter x at site y as specified in neuron{n}.record
                                    end
                                    
                                    fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                    if params.cvode
                                        fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                    end
                                end
                            end
                            neuron{n}.record{t}{r,4} = realrecs;
                        case 'pp'
                            if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                for in =  1:size(realrecs,1)%1:numel(neuron{n}.record{t}{r,1})
                                    fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                    fprintf(ofile,sprintf('rec.label("%s of %s Point Process #%d at location %06.4f of section %d of cell %d")\n', neuron{n}.record{t}{r,2} , neuron{n}.pp{t}{neuron{n}.record{t}{r,1},1}, neuron{n}.record{t}{r,1}, minterf{t}(find(minterf{t}(:,1) == neuron{n}.pp{t}{neuron{n}.record{t}{r,1},2},1,'first'),[4 2]) ,t-1) ); % label the vector for plotting
                                    if params.cvode
                                        fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&ppList.o(%d).%s,rec,rect)}\n',t-1, realrecs(in,1),sum(ppnum(1:t-1))+in-1, neuron{n}.record{t}{r,2} ) ); % record the parameter x at site y as specified in neuron{n}.record
                                    else
                                        fprintf(ofile,sprintf('io = rec.record(&ppList.o(%d).%s,tvec)\n',sum(ppnum(1:t-1))+in-1, neuron{n}.record{t}{r,2} ) ); % record the parameter x at site y as specified in neuron{n}.record
                                    end
                                    fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                    if params.cvode
                                        fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                    end
                                end
                            end
                        case 'stim'
                            if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                for in = 1: 1:size(realrecs,1)%numel(neuron{n}.record{t}{r,1})
                                    fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                    fprintf(ofile,sprintf('rec.label("%s of %s electrode %d at location %06.4f of section %d of cell %d")\n', neuron{n}.record{t}{r,2} , neuron{n}.stim{t}{neuron{n}.record{t}{r,1},1}, neuron{n}.record{t}{r,1}, minterf{t}(find(minterf{t}(:,1) == neuron{n}.stim{t}{neuron{n}.record{t}{r,1},2},1,'first'),[4 2]) ,t-1) ); % label the vector for plotting
                                    if params.cvode
                                        fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec {io = cvode.record(&stimList.o(%d).%s,rec,rect)}\n',t-1, realrecs(in,1),sum(stimnum(1:t-1))+in-1, neuron{n}.record{t}{r,2} ) ); % record the parameter x at site y as specified in neuron{n}.record, CAUTION! Absolutely necessary to temporally access section! Cvode needs this!!!
                                    else
                                        fprintf(ofile,sprintf('io = rec.record(&stimList.o(%d).%s,tvec)\n',sum(stimnum(1:t-1))+in-1, neuron{n}.record{t}{r,2} ) ); % record the parameter x at site y as specified in neuron{n}.record
                                    end
                                    
                                    fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                    if params.cvode
                                        fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                    end
                                end
                            end
                        case 'artificial'
                            if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                fprintf(ofile,sprintf('rec = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                fprintf(ofile,sprintf('rec.label("%s of artificial cell %s (cell #%d)")\n', neuron{n}.record{t}{r,2} , tree{t}.artificial, t-1) ); % label the vector for plotting
%                                 if strcmpi(neuron{n}.record{t}{r,2},'on')
%                                     fprintf(ofile,sprintf('nilcon = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',t-1,0.5) );    % for art. cells, make netcon with threshold 0.5
%                                     fprintf(ofile,sprintf('io = rec.record(&nilcon)\n')fprintf(ofile,'io = APC.record(APCrecList.o(APCrecList.count()-1))\n');
%                                 else
                                    if params.cvode
                                        fprintf(ofile,sprintf('rect = new Vector(%f)\n',(params.tstop-params.tstart)/params.dt+1 ) );    % create new recording vector
                                        fprintf(ofile,sprintf('io = cvode.record(&cellList.o(%d).cell.%s,rec,rect)\n',t-1, neuron{n}.record{t}{r,2} ) );  % record the parameter x of artificial cell t-1
                                    else
                                        fprintf(ofile,sprintf('io = rec.record(&cellList.o(%d).cell.%s,tvec)\n', t-1, neuron{n}.record{t}{r,2} ) ); % record the parameter x of artificial cell t-1
                                    end
%                                 end
                                fprintf(ofile,'io = recList.append(rec)\n\n' );  %append recording vector to recList
                                if params.cvode
                                    fprintf(ofile,'io = rectList.append(rect)\n\n' );  %append time recording vector to recList
                                end
                            end
                    end
                    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                        fprintf(ofile,'\n');
                    end
                end
            end
        end
        if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
            fprintf(ofile, 'objref rec\n');
            fprintf(ofile, 'objref rect\n');
        end
    end
    
    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
        %!%! there might be problems with cvode (not adjusted yet)
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define APCount sites *****\n');
        if isfield(neuron{n},'APCount')
            for t = 1:numel(tree)
                if numel(neuron{n}.APCount) >= t && ~isempty(neuron{n}.APCount{t})   % if a recording site was defined for  this tree
                    for r = 1: size(neuron{n}.APCount{t},1)
                        if ~isfield(tree{t},'artificial')
                            inode = find(minterf{t}(:,1) == neuron{n}.APCount{t}(r,1),1,'first');    %find the index of the node in minterf
                            fprintf(ofile,sprintf('cellList.o(%d).allregobj.o(%d).sec',t-1,minterf{t}(inode,2) ) );    % corresponding section of node
                            fprintf(ofile,sprintf('{APC = new APCount(%f)\n',minterf{t}(inode,3) ) );    % make APCCount at position x
                            fprintf(ofile,sprintf('APC.thresh = %f\n',neuron{n}.APCount{t}(r,2) ) ); % set threshold of APCount [mV]
                            %             if neuron{n}.APCount{t}{r,3}
                        else
                            fprintf(ofile,sprintf('APC = new NetCon(cellList.o(%d).cell,nil,%g,0,5)\n',t-1,neuron{n}.APCount{t}(r,2) ) );    % for art. cells, make netcon with threshold
                        end
                        fprintf(ofile,'APCrec = new Vector()\n');
                        fprintf(ofile,'io = APCrecList.append(APCrec)\n');
                        fprintf(ofile,'io = APC.record(APCrecList.o(APCrecList.count()-1))\n');
                        %             end
                        
                        if ~isfield(tree{t},'artificial')
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
        end
        fclose(ofile);
    end
    %% write init_play.hoc
    if params.changed.play || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'init_play.hoc') ,'wt');   %open record hoc file in write modus
        fprintf(ofile,'\n\n');
        fprintf(ofile,'// ***** Define play sites *****\n');
        if isfield(neuron{n},'play')
            for t = 1:numel(tree)
                if numel(neuron{n}.play) >= t &&~isempty(neuron{n}.play{t})  && ~isfield(tree{t},'artificial')   % if a playing site was defined for  this tree
                    for p = 1: size(neuron{n}.play{t},1)
                        inode = find(minterf{t}(:,1) == neuron{n}.play{t}{p,1},1,'first');    %find the index of the node in minterf
                        fprintf(ofile,sprintf('playt = new Vector(%f)\n',length(neuron{n}.play{t}{p,3}) ) );    % create new playing time vector
                        %a file needs to be created to temporally save the vector so
                        %NEURON can read it in. otherwise it would be necessary to
                        %print the whole vector into the hoc file. alternatively i
                        %could give a file name where the vector lies so it is not
                        %written each time cn is called...
                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('plt_%s_at_%d_cell_%d.dat', neuron{n}.play{t}{p,2} , inode ,t-1)),'w');
                        fprintf(f,'%g ', neuron{n}.play{t}{p,3}(1:end-1));
                        fprintf(f,'%g\n', neuron{n}.play{t}{p,3}(end));
                        fclose(f);
                        fprintf(ofile,'f = new File()');
                        fprintf(ofile,sprintf('f.ropen("plt_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}{p,2} , inode ,t-1));  %vector file is opened
                        fprintf(ofile,'playt.scanf(f)');    % file is read into time vector
                        fprintf(ofile,'io = f.close()');     %file is closed
                        fprintf(ofile,'io = playtList.append(playt)\n\n' );  %append playing time vector to playtList
                        
                        fprintf(ofile,sprintf('play = new Vector(%f)\n',length(neuron{n}.play{t}{p,4}) ) );    % create new playing vector
                        f = fopen(fullfile(exchfolder,thisfolder,sprintf('pl_%s_at_%d_cell_%d.dat', neuron{n}.play{t}{p,2} , inode ,t-1)),'w');
                        fprintf(f,'%g ', neuron{n}.play{t}{p,4}(1:end-1));
                        fprintf(f,'%g\n', neuron{n}.play{t}{p,4}(end));
                        fclose(f);
                        fprintf(ofile,'f = new File()');
                        fprintf(ofile,sprintf('f.ropen("pl_%s_at_%d_cell_%d.dat")\n', neuron{n}.play{t}{p,2} , inode ,t-1));  %vector file is opened
                        fprintf(ofile,'play.scanf(f)');     % file is read into play vector
                        fprintf(ofile,'io = f.close()');   %file is closed
                        fprintf(ofile,sprintf('play.label("playing %s at node %d of cell %d")\n', neuron{n}.play{t}{p,2} , inode ,t-1) ); % label the vector for plotting
                        fprintf(ofile,sprintf('play.play(&cellList.o(%d).allregobj.o(%d).sec.%s(%f),playtList.o(playtList.count()-1),%d)\n',t-1,minterf{t}(inode,2), neuron{n}.play{t}{p,2}, minterf{t}(inode,3), neuron{n}.play{t}{p,5} ) ); % play the parameter x at site y as specified in neuron{n}.play
                        fprintf(ofile,'io = playList.append(play)\n\n' );  %append playing vector to playList
                        
                    end
                    fprintf(ofile,'\n');
                end
            end
            fprintf(ofile, 'objref playt\n');
            fprintf(ofile, 'objref play\n');
        end
        fclose(ofile);
    end
    
    %% write save_rec.hoc
    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
        ofile = fopen(fullfile(exchfolder,thisfolder,'save_rec.hoc') ,'wt');   %open record hoc file in write modus
        fprintf(ofile,'// * Write Recordings to Files *\n');
    end
    if isfield(neuron{n},'record')
        out{n}.record = cell(1,numcell);   % initialize output of cn
        
        c = 0;
        
        for t = 1:numel(tree)
            if numel(neuron{n}.record) >= t && ~isempty(neuron{n}.record{t})
                for r = 1: size(neuron{n}.record{t},1)
                    if isfield(tree{t},'artificial')
                        rectype = 'artificial';
                        if size(neuron{n}.record{t},2) < 2 || isempty(neuron{n}.record{t}{r,2})
                            neuron{n}.record{t}{r,2} = neuron{n}.record{t}{r,1};                    % if parameter to record was defined in the first entry
                        end
                    elseif size(neuron{n}.record{t},2)<3 || isempty(neuron{n}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                        rectype = 'node';
                    else
                        rectype = neuron{n}.record{t}{r,3};
                        if any(strcmpi(rectype,{'PP','Syn'}))        % in this case also the electrodes can be put to the Point Processes
                            rectype = 'pp';
                        elseif strcmpi(rectype,{'Stim'})
                            rectype = 'stim';
                        end
                    end
                    switch rectype
                        case 'node'
                            for in = 1:size(neuron{n}.record{t}{r,4},1)
                                fname = sprintf('cell%d_sec%d_loc%06.4f_%s',t-1, neuron{n}.record{t}{r,4}(in,1), neuron{n}.record{t}{r,4}(in,2), neuron{n}.record{t}{r,2} );
                                if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                    fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                    fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname) )  );  % open file for this vector with write perm.
                                    fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                    fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                    if params.cvode
                                        fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                        fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname) )  );  % open file for this vector with write perm.
                                        fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                        fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                    end
                                end
                                c= c+1;
                                noutfiles = noutfiles +1;
                                %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{n}.record{t}{r,1},neuron{n}.record{t}{r,2} ) , 'record' ,  t , neuron{n}.record{t}{r,2} ,neuron{n}.record{t}{r,1} };
                                readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{n}.record{t}{r,2} ,neuron{n}.record{t}{r,4}(in,1) ,neuron{n}.record{t}{r,4}(in,2) };
                            end
                        case {'pp','stim'}
                           for in = 1:numel(neuron{n}.record{t}{r,1})
                               if strcmp(rectype,'pp')
                                   fname = sprintf('cell%d_%s_%d_%s',t-1, neuron{n}.pp{t}{neuron{n}.record{t}{r,1},1}, neuron{n}.record{t}{r,1}(in) , neuron{n}.record{t}{r,2} );
                               else
                                   fname = sprintf('cell%d_%s_%d_%s',t-1, neuron{n}.stim{t}{neuron{n}.record{t}{r,1},1}, neuron{n}.record{t}{r,1}(in) , neuron{n}.record{t}{r,2} );
                               end
                               if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                   fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                   fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname))  );  % open file for this vector with write perm.
                                   fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                   fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                   if params.cvode
                                       fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                       fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname))  );  % open file for this vector with write perm.
                                       fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                       fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                   end
                               end
                                c= c+1;
                                noutfiles = noutfiles +1;
                                %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{n}.record{t}{r,1},neuron{n}.record{t}{r,2} ) , 'record' ,  t , neuron{n}.record{t}{r,2} ,neuron{n}.record{t}{r,1} };
                                readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{n}.record{t}{r,2} };
                           end
                        case 'artificial'
                            fname = sprintf('cell%d_%s',t-1, neuron{n}.record{t}{r,2} );
                            if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                                fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s.dat',fname))  );  % open file for this vector with write perm.
                                fprintf(ofile,sprintf('io = recList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                if params.cvode
                                    fprintf(ofile,'f = new File()\n');      %create a new filehandle
                                    fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,sprintf('%s_tvec.dat',fname))  );  % open file for this vector with write perm.
                                    fprintf(ofile,sprintf('io = rectList.o(%d).printf(f, "%%%%-20.20g\\\\n")\n', c ) );    % print the data of the vector into the file
                                    fprintf(ofile,'io = f.close()\n');   %close the filehandle
                                end
                            end
                            c= c+1;
                            noutfiles = noutfiles +1;
                            %                     readfiles{noutfiles} = {sprintf('cell%d_node%d_%s.dat',t-1,neuron{n}.record{t}{r,1},neuron{n}.record{t}{r,2} ) , 'record' ,  t , neuron{n}.record{t}{r,2} ,neuron{n}.record{t}{r,1} };
                            readfiles{noutfiles} = {n, sprintf('%s.dat',fname) , rectype ,  t , neuron{n}.record{t}{r,2} };
                    end
                end
            end
        end
        if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
            fprintf(ofile,'\n');
        end
    end
    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
        fprintf(ofile,'// * Write APCounts to Files *\n');
    end
    if isfield(neuron{n},'APCount')
        out{n}.APCtimes = cell(1,numcell);   % initialize output of cn
        
        c=0;
        for t = 1:numel(tree)
            if numel(neuron{n}.APCount) >= t && ~isempty(neuron{n}.APCount{t})     % if a recording site was defined for  this tree
                for r = 1: size(neuron{n}.APCount{t},1)
                    fname = sprintf('cell%d_node%d_APCtimes.dat',t-1,neuron{n}.APCount{t}(r,1) );
                    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                        fprintf(ofile,'f = new File()\n');      %create a new filehandle
                        fprintf(ofile,sprintf('io = f.wopen("%s//%s//%s")\n',nrn_exchfolder,thisfolder,fname) );  % open file for this vector with write perm.
                        fprintf(ofile,sprintf('io = APCrecList.o(%d).printf(f, "%%%%-20.10g")\n', c ) );    % print the data of the vector into the file
                        fprintf(ofile,'io = f.close()\n');   %close the filehandle
                    end
                    c= c+1;
                    noutfiles = noutfiles +1;
                    readfiles{noutfiles} = {n, fname , 'APCtimes' , t , neuron{n}.APCount{t}(r,1) };
                    
                end
                if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
                    fprintf(ofile,'\n');
                end
            end
        end
    end
    if params.changed.rec || params.changed.morph     %rewrite only if something has changed influencing this file
        fclose(ofile);
    end
    
    if strfind(options,'-cl') %transfer files to server
        filenames = {interf_file,'init_cells.hoc','init_mech.hoc','init_pp.hoc','init_con.hoc','init_stim.hoc','init_rec.hoc','save_rec.hoc','init_play.hoc'}; %'init_pas.hoc',
        localfilename = cell(max(sum(structfun(@(x) x,params.changed))-1+params.changed.rec,params.changed.morph*10),1);
        remotefilename = cell(max(sum(structfun(@(x) x,params.changed))-1+params.changed.rec,params.changed.morph*10),1);
        m = 1;
        if params.changed.basic || params.changed.lib || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{1});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{1});
            m = m + 1;
        end
        if params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{2});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{2});
            m = m + 1;
        end
        if  params.changed.mech || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{4});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{4});
            m = m + 1;
        end
        if params.changed.pp || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{5});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{5});
            m = m + 1;
        end
        if params.changed.con || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{6});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{6});
            m = m + 1;
        end
        if params.changed.stim || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{7});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{7});
            m = m + 1;
        end
        if params.changed.rec || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{8});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{8});
            m = m + 1;
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{9});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{9});
            m = m + 1;
        end
        if params.changed.play || params.changed.morph
            localfilename{m} = fullfile(exchfolder,thisfolder,filenames{10});
            remotefilename{m} = sprintf('%s/%s/%s',nrn_exchfolder,thisfolder,filenames{10});
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
        fprintf(sprintf('HOC writing time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
end

%% Execute NEURON
if ~isempty(strfind(options,'-cl'))
    warndlg('Multi-Thread and server execution not yet implemented')
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
            % NOT CHECKED YET
            s = find(simids==1);
            for ss = 1:numel(s)
%                 if ~isempty(result{s(ss)}) && ~isempty(strfind(result{s(ss)},'error')) || ~isempty(strfind(result{s(ss)},'near line'))
%                     ret = regexp(result{s(ss)},'\n');
%                     er =  strfind(result{s(ss)},'error');
%                     if isempty(er)
%                         er = strfind(result{s(ss)},'near line');
%                         er = er(1);
%                         ind = find(ret < er,1,'last')-1;
%                     else
%                         er = er(1);
%                         ind = find(ret < er,1,'last');
%                     end
%                     result{s(ss)} = result{s(ss)}(ret(ind)+1:end);
%                     errordlg(sprintf('An error occurred during NEURON execution:\n******************************\n%s\n******************************\nDue to that, m2n does not return an output!',result{s(ss)}))
%                     r = find(simids==0,1,'first');  % find next not runned simid
%                     if ~isempty(r)
%                         [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
%                         simids(r) = 1;          % mark this as running
%                         simids(s(ss)) = 3;
%                     end
% 
%                 end
            
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
                        errordlg(sprintf('There was an error in Simulation %d:\n******************************\n%s\n******************************\nDue to that m2n has no output to that Simulation.',s(ss),txt));
                        simids(s(ss)) = 3;
                        r = find(simids==0,1,'first');  % find next not runned simid
                        if ~isempty(r)
                            [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                            simids(r) = 1;          % mark this as running
                        end
                        %                 elseif exist(fullfile(exchfolder,sprintf('sim%d',s(ss)),'NeuronLogFile.txt'),'file') == 2
                        %                     f = fopen(fullfile(exchfolder,sprintf('sim%d',s(ss)),'NeuronLogFile.txt'));
                        %                     txt = textscan(f,'%s');
                        % %                     txt = fscanf(f,'%c');
                        %                     fclose(f);
                        %                     for x = 1:numel(txt)
                        %                         if strfind(txt{x},{'fail','error'})
                        %                             'g'
                        %                             r = find(simids==0,1,'first');  % find next not runned simid
                        %                             if ~isempty(r)
                        %                                 [jobid(r),tim] = exec_neuron(r,exchfolder,nrn_exchfolder,interf_file,params,options);          % start new simulation
                        %                                 simids(r) = 1;          % mark this as running
                        %                             end
                        %                         end
                        %                     end
                    end
                end
            end
            pause(0.1);
        end
        if ~isempty(strfind(options,'-w'))
            if ishandle(w)
                waitbar(sum(simids>1)/numel(simids),w)
            else
                errordlg('Waitbar was closed, m2n stopped continuing. If accidently, retry.')
                out.error = 1;
                return
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
    %
    %         dos('exit');  % exit NEURON if it was defined so in the parameters
    %     end
    if ~isempty(strfind(options,'-d'))
        tim = toc(timm);
        fprintf(sprintf('NEURON execute time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
    if isempty(strfind(options,'-q'))
        display('NEURON finished... loading data...')
    end
    if strfind(options,'-d')
        tim = tic;
    end
    if strfind(options,'-cl')
        outputnames = cellfun(@(x) strcat(nrn_exchfolder,sprintf('/sim%d/',x{1}),x{2}),readfiles,'UniformOutput',0);  % extract filenames
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
    for f = 1:noutfiles
        if simids(readfiles{f}{1}) == 2    % if there was no error during simulation
            fn = fullfile(exchfolder,sprintf('sim%d',readfiles{f}{1}),readfiles{f}{2});
            switch readfiles{f}{3}
                case {'node' , 'artificial' , 'pp', 'stim'}
                    readfiles{f}{8} = load(fn,'-ascii');    %temporary loading of file. association is done below
                    if params.cvode
                        if params.use_local_dt  % if yes, dt was different for each cell, so there is more than one time vector
                            if numel(out{readfiles{f}{1}}.t) < readfiles{f}{4} || isempty(out{readfiles{f}{1}}.t{readfiles{f}{4}})
                                out{readfiles{f}{1}}.t{readfiles{f}{4}} = load(strcat(fn(1:end-4),'_tvec',fn(end-3:end)),'-ascii');    %loading of one time vector file per cell (sufficient)
                            end
                        elseif ~ isfield(out{readfiles{f}{1}},'t')       % if it has been loaded in a previous loop
                            out{readfiles{f}{1}}.t = load(strcat(fn(1:end-4),'_tvec',fn(end-3:end)),'-ascii');    %loading of one time vector file at all (sufficient)
                        end
                        out{readfiles{f}{1}}.t(find(diff(out{readfiles{f}{1}}.t,1) == 0) + 1) = out{readfiles{f}{1}}.t(find(diff(out{readfiles{f}{1}}.t,1) == 0) + 1) + 1e-10;  % add tiny time step to tvec to avoid problems with step functions
                        %   not necessary yet       readfiles{f}{9} = load(strcat(fn(1:end-4),'_tvec',fn(end-3:end)),'-ascii');    %temporary loading of file. association is done below
                    end
                case 'APCtimes'
                    out{readfiles{f}{1}}.(readfiles{f}{3}){readfiles{f}{4}}{readfiles{f}{5}} = load(fn,'-ascii');
                otherwise
                    errordlg(sprintf('Data "%s" not specified for output',readfiles{f}{2}))
            end
        else
            out{readfiles{f}{1}}.error = 1;
        end
    end
    if isfield(neuron{n},'record')
        for n = 1:numel(neuron)
            if simids(n) == 2       % if there was no error during simulation
                for t = 1:numel(tree)
                    if numel(neuron{n}.record) >= t && ~isempty(neuron{n}.record{t})  % if a recording site was defined for  this tree
                        for r = 1: size(neuron{n}.record{t},1)     %go through all set recordings
                            if isfield(tree{t},'artificial')
                                rectype = 'artificial';
                                if size(neuron{n}.record{t},2) < 2 || isempty(neuron{n}.record{t}{r,2})
                                    neuron{n}.record{t}{r,2} = neuron{n}.record{t}{r,1};                    % if parameter to record was defined in the first entry
                                end
                            elseif size(neuron{n}.record{t},2)<3 || isempty(neuron{n}.record{t}{r,3})       % check if recording should be a parameter in a section, or a point process (divided in pps and electrodes)
                                rectype = 'node';
                            else
                                rectype = neuron{n}.record{t}{r,3};
                                if any(strcmp(rectype,{'PP','Pp','syn','Syn','stim','Stim','STIM'}))
                                    rectype = 'pp';
                                end
                            end
                            
                            switch rectype
                                case 'node'
                                    for in = 1:numel(neuron{n}.record{t}{r,1})  %go through all nodes in this defined recording
                                        inode = find(minterf{t}(:,1) == neuron{n}.record{t}{r,1}(in),1,'first');    %find the index of the node in minterf
                                        %the correct file for this node is searched in the
                                        %temporally loaded files and copied
                                        onlyrecords = cellfun(@(x) strcmpi(x{3},'node'),readfiles);   % don't know any more why that complicated but never change a running system^^
                                        thisfile = find(cumsum(onlyrecords) == find(cellfun(@(x) n==x{1} & strcmpi(x{5},neuron{n}.record{t}{r,2}) & x{4} == t & x{6} == minterf{t}(inode,2) & x{7} == minterf{t}(inode,4) ,readfiles(onlyrecords))),1,'first');
                                        out{n}.record{t}.(neuron{n}.record{t}{r,2}){neuron{n}.record{t}{r,1}(in)} = readfiles{thisfile}{8};
                                        % alternatively only give pointer and give temp
                                        % files as output, too...
                                    end
                                case 'pp'
                                    for in = 1:numel(neuron{n}.record{t}{r,1})
                                        thisfile = find(cellfun(@(x) n==x{1} & any(strcmpi(x{3},{'pp','stim'})) & strcmpi(x{5},neuron{n}.record{t}{r,2}) & x{4} == t ,readfiles),1,'first');
                                        out{n}.record{t}.(neuron{n}.record{t}{r,2}){neuron{n}.record{t}{r,1}(in)} = readfiles{thisfile}{8};
                                    end
                                case 'artificial'
                                    thisfile = find(cellfun(@(x) n==x{1} & strcmpi(x{3},'artificial') & strcmpi(x{5},neuron{n}.record{t}{r,2}) & x{4} == t ,readfiles),1,'first');
                                    out{n}.record{t}.(neuron{n}.record{t}{r,2}) = readfiles{thisfile}{8};
                            end
                            
                            
                        end
                    end
                end
            end
        end
    end
    if ~isempty(strfind(options,'-w')) 
        close(w)
    end
    if isempty(strfind(options,'-q'))
        display('data sucessfully loaded')
    end
    if strfind(options,'-d')
        tim = toc(tim);
        fprintf(sprintf('Data loading time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))
    end
else
    out = [];
end

if nargout < 4 && orderchanged
        warndlg('Caution, the node order of some trees had to be changed! Sort your trees with "sort_tree" to obtain the correct results','Node order change!')
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
if ~isempty(strfind(options,'-cl')) && numel(answer) == 1 && isempty(strfind(options,'-f'))
    str = regexp(answer{1},'[0-9]*','match');
    ncount = cellfun(@numel,str);
    [nix, ind] = max(ncount);
    jobid = str2double(str{ind});
    % there might be error if many jobs are run, because answer might not
    % be 1
else
    jobid = NaN;
end
end

function minterf = make_nseg(tree,minterf,params,mech)
%does the same as the d_lambda procedure in NEURON
%necessary to find nearest segment which will be calculated
if ischar(params.nseg) && strcmpi(params.nseg,'dlambda')
    dodlambda = 1;
    pl = PL_tree(tree);     % path length of tree..%%%%not anymore %%%%add zero because neuron_template_tree adds one tiny segment at root
    D =  tree.D; %%%not anymore%%%%%add zero because neuron_template_tree adds one tiny segment at root
    freq = 100;
    %     if ~isempty(mech) && all(~isnan(mech))
    %         Ra = pas(2);
    %         cm = pas(1);
    %     elseW
    %
    %     end
    if isfield(params,'d_lambda')
        d_lambda = params.d_lambda;
    else
        
        d_lambda = 0.1;
    end
else
    dodlambda = 0;
end

for sec = 0:max(minterf(:,2))  %go through all sections
    
    secstart = find(minterf(:,2) == sec & minterf(:,3) == 0);
    secend = find(minterf(:,2) == sec & minterf(:,3) == 1);
    if dodlambda
        secnodestart = minterf(secstart,1);
        if secnodestart == 0  % this means the current section is the tiny section added for neuron... this should have nseg = 1
            secnodestart = 1;
            flag = true;
        else
            flag = false;
        end
        secnodestart2 = minterf(secstart+1,1);
        if isfield(tree,'rnames') && isfield(mech,tree.rnames{tree.R(secnodestart)}) && isfield(mech.(tree.rnames{tree.R(secnodestart)}),'pas') && all(isfield(mech.(tree.rnames{tree.R(secnodestart)}).pas,{'Ra','cm'}))
            Ra = mech.(tree.rnames{tree.R(secnodestart)}).pas.Ra;
            cm = mech.(tree.rnames{tree.R(secnodestart)}).pas.cm;
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
        %         a=zeros(numel(tree.X),1)
        %         a(secnodestart:secnodeend)=1
        %         plot_tree(tree,a)
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
        %         fprintf('%g\n',(L/(d_lambda*lambda_f)+0.9))
    else
        nseg = params.nseg;
    end
    %     fprintf('%d\n',nseg);
    if isfield(params,'accuracy') && params.accuracy == 1 %triple nseg if accuracy is necessary
        nseg = 3 * nseg;
    end
    %     fprintf('%d\n',nseg)
    pos = (2 * (1:nseg) - 1) / (2*nseg);    %calculate positions
    for in = secstart:secend
        [nix,ind] = min(abs(minterf(in,3) - pos));   %find position to node which is closest to next segment location
        minterf(in,4) = pos(ind);                % hier evt ausnahme für anfang und ende der section (=0)
    end
end

end