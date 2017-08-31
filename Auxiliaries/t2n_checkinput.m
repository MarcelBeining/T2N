function [neuron,tree,usestreesof,nocell,exchfolder] = t2n_checkinput(neuron,tree,options)
% This function checks the neuron structure for correct definition of the
% used morphologies and returns info about it

% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
% options           (optional) option string of t2n
%
% OUTPUTS
% tree              corrected tree cell array
% neuron            corrected neuron structure
% usestreesof       points to the neuron entry/instance from which the tree
%                   definitions are taken from
% nocell            Boolean if neuron input was a structure or cell array
% exchfolder        name for exchfolder that was possibly found in the
%                   neuron structure
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 4 || isempty(options)
    options = '';
end
if ~exist(fullfile(pwd,'morphos','hocs'),'dir')
    mkdir(fullfile(pwd,'morphos','hocs'));
end

%% check neuron
nocell = false;
if isstruct(neuron)         % transform structure neuron to cell neuron
    if numel(neuron) == 1
        neuron = {neuron};
        nocell = true;
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
if  ~all(cellfun(@(x) isfield(x,'NID'),tree))
    doit = 1;
else
    NIDs = unique(cellfun(@(x) x.NID,tree,'UniformOutput',0));  % improves speed if many same cells are used
    doit = 0;
end
if doit || ~all(cellfun(@(x) exist(fullfile(pwd,'morphos','hocs',strcat(x,'.hoc')),'file'),NIDs))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet or hoc file is missing! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        ind = cellfun(@(x) ~isfield(x,'NID'),tree) ;
        if ~all(ind)
            ind(~ind) = ~cellfun(@(x) exist(fullfile(pwd,'morphos','hocs',strcat(x.NID,'.hoc')),'file'),tree(~ind));
        end
        tree(ind) = t2n_writeTrees(tree(ind));
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
            else
                n = find(bool);
                flag = true;
            end
        end
    case 0      % if no trees are given, use trees that are given to t2n in their order...
        thesetrees = repmat({1:numel(tree)},numel(neuron),1);
        usestreesof = ones(numel(neuron),1);
    case numel(neuron)
        for n = 1:numel(neuron)
            x = t2n_getref(n,neuron,'tree');
            if ~isnan(x)
                thesetrees{n} = unique(neuron{x}.tree);
                usestreesof(n) = x;
            elseif isempty(x)
                thesetrees{n} = 1:numel(tree);%repmat({1:numel(tree)},numel(neuron),1);
                usestreesof(n) = 1;%ones(numel(neuron),1);
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

if isfield(neuron{1}.params,'exchfolder')
    exchfolder = neuron{1}.params.exchfolder;
else
    exchfolder = [];
end
for n = 1:numel(neuron)
    neuron{n}.tree = thesetrees{n};
    
    %% check for several standard parameter and initialize default value if not set
    refP = t2n_getref(n,neuron,'params');
    if n == refP % only check if current instance has its own parameter struct
        if ~isfield(neuron{n}.params,'parallel')
            neuron{n}.params.parallel = false;
        end
        if ~isfield(neuron{n}.params,'cvode')
            neuron{n}.params.cvode = false;
        end
        if ~isfield(neuron{n}.params,'use_local_dt')
            neuron{n}.params.use_local_dt = 0;
        end
        if ~isfield(neuron{n}.params,'nseg')
            neuron{n}.params.nseg = 'dlambda';
            disp('Number of segments or nseg rule not set in neuron{n}.params.nseg. Dlambda rule will be applied')
        elseif strcmpi(neuron{n}.params.nseg,'d_lambda')
            neuron{n}.params.nseg = 'dlambda';
        end
        if ~isfield(neuron{n}.params,'tstart')
            neuron{n}.params.tstart = 0;
        end
        if ~isfield(neuron{n}.params,'tstop')
            neuron{n}.params.tstop = 200;
            disp('Tstop not defined in neuron{n}.params.tstop. Default value of 200 ms is applied.')
        end
        if ~isfield(neuron{n}.params,'dt')
            neuron{n}.params.dt = 0.025;
            if ~neuron{n}.params.cvode
                disp('Time step not defined in neuron{n}.params.dt. Default value of 0.025 ms is applied.')
            end
        end
        if ~isfield(neuron{n}.params,'accuracy')
            neuron{n}.params.accuracy = 0;
        end
        if ~isfield(neuron{n}.params,'skiprun')
            neuron{n}.params.skiprun = false;
        end
        if ~isfield(neuron{n}.params,'q10')
            neuron{n}.params.q10 = false;
        end
        if ~isfield(neuron{n}.params,'prerun')
            neuron{n}.params.prerun = false;
        end
        if neuron{n}.params.cvode && isnumeric(neuron{n}.params.dt) && (~isnan(neuron{n}.params.dt) || ~isempty(neuron{n}.params.dt))
            warning ('t2n:cvode', 'Dt is set but cvode is active. Dt will be ignored');
        end
    end

end
end