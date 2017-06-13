function [tree,params,neuron,thesetrees,usestreesof] = t2n_checkinput(tree,params,neuron)
% This function checks the neuron structure for correct definition of the
% used morphologies and returns info about it

% INPUTS
% tree              tree cell array with morphologies (see documentation)
% params            t2n parameter structure (see documentation)
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
%
% OUTPUTS
% tree              corrected tree cell array
% params            corrected parameter structure
% neuron            corrected neuron structure
% thesetrees        returns for each cell array entry in neuron the index
%                   to the trees that are used in this neuron instance
% usestreesof       points to the neuron entry/instance from which the tree
%                   definitions are taken from
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

%% check for several standard parameter and initialize default value if not set
if ~isempty(params)
    if ~isfield(params,'parallel')
        params.parallel = false;
    end
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
    if ~isfield(params,'path')
        params.path = regexprep(pwd,'\\','/');
        display('No standard path was given. Current folder is used instead');
    else
        params.path = regexprep(params.path,'\\','/');
    end
    if strcmpi(params.path(end),'\')
        params.path = params.path(1:end-1);
    end
    if params.cvode && isnumeric(params.dt)
        warning ('t2n:cvode', 'Dt is set but cvode is active. Dt will be ignored');
    end
    if ~isfield(params,'neuronpath')
        if ispc
            warning('Path to neuron not given in params.neuronpath! Default path chosen (which might be wrong)')
        end
        params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';  % default path to neuron exe
    end
    if ~isfield(params,'morphfolder') % check for morphology folder
        if exist(fullfile(nrn_path,'morphos'),'dir')
            params.morphfolder = 'morphos';
        else
            error('Please give morphology folder in params.morphfolder or create standard morphology folder "morphos"');
        end
    end
    parflag = true;
else
    parflag = false;
end
%% check neuron
params.nocell = false;
if isstruct(neuron)         % transform structure neuron to cell neuron
    if numel(neuron) == 1
        neuron = {neuron};
        params.nocell = true;
    else
        neuron = arrayfun(@(y) {y},neuron);
    end
end

%% check tree
if iscell(tree) && iscell(tree{1})
    tree = tree{1};
elseif isstruct(tree)
    tree = {tree};
end
for t = 1:numel(tree)
    if isfield(tree{t},'artificial') && ~isfield(tree{t},'NID')
        tree{t}.NID = strcat('cell_',tree{t}.artificial);           % artificial cells only need one morph hoc file which is named cell_ + the name of the artificial cell..
    end
end
NIDs = unique(cellfun(@(x) x.NID,tree,'UniformOutput',0));  % improves speed if many same cells are used
if parflag && (~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(params.path,params.morphfolder,strcat(x,'.hoc')),'file'),NIDs)))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet or hoc file is missing! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        ind = cellfun(@(x) ~isfield(x,'NID'),tree) ;
        if ~all(ind)
            ind(~ind) = ~cellfun(@(x) exist(fullfile(params.path,params.morphfolder,strcat(x.NID,'.hoc')),'file'),tree(~ind));
        end
        
        tree(ind) = t2n_writeTrees(tree(ind),params,'',options);
    else
        error('Aborted');
    end
end

%% check tree/neuron consistency
thesetrees = cell(numel(neuron),1);
usestreesof = zeros(numel(neuron),1);
flag = false;
bool = cellfun(@(y) isfield(y,'tree'),neuron);
if all(cellfun(@(y) ischar(y.tree),neuron(bool))) % no trees defined (only due to t2n_as function). make bool = 0
    bool(:) = 0;
end
nempty = cellfun(@isempty,neuron);
if any(nempty)
    error('The defined simulation #%d is empty, please check\n',find(nempty))
end
switch sum(bool)
    case 1   % use that treeids defined in that one simulation
        if isnumeric(neuron{bool}.tree)
            thesetrees = repmat({unique(neuron{bool}.tree)},numel(neuron),1);
            usestreesof = repmat(find(bool),numel(neuron),1);
        else
            x = t2n_getref(1,neuron,'tree');
            if isempty(x) % means it refers to itself (maybe due to usage of t2n_as)..use normal trees..
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
    case 0      % if no trees are given, use trees that are given to t2n in their order...
        thesetrees = repmat({1:numel(tree)},numel(neuron),1);
        usestreesof = ones(numel(neuron),1);
        for n = 1:numel(neuron)
            neuron{n}.tree = thesetrees{n};
        end
    case numel(neuron)
        for n = 1:numel(neuron)
            x = t2n_getref(n,neuron,'tree');
            if ~isnan(x)
                thesetrees{n} = unique(neuron{x}.tree);
                usestreesof(n) = x;
            elseif isempty(x)
                thesetrees{n} = 1:numel(tree);%repmat({1:numel(tree)},numel(neuron),1);
                usestreesof(n) = 1;%ones(numel(neuron),1);
                %                 for n = 1:numel(neuron)
                neuron{n}.tree = thesetrees{n};
                %                 end
            else
                flag = true;
                break
            end
        end
    otherwise  % if more than one are given, t2n cannot know which trees you want
        n = find(bool);
        flag = true;
end
if flag
    error('Error in neuron{%d}.tree, please check\n',n)
end
if ~parflag
    params = [];
end
end