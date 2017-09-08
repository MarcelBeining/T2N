function  [GIDs,neuron,mindelay] = t2n_getGIDs(neuron,tree,thesetrees)
% This function prepares the neuron structure to be used for the parallel
% NEURON environment. Also it returns the global IDs, which are necessary
% for the parallel environment. The parallel NEURON feature is not yet
% implemented completely!
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
% thesetrees        (optional) if not all trees of 'tree' are used, this is
%                   the index array to the used ones
%
% OUTPUTS
% GIDs              global IDs of each tree for initializing NEURON in
%                   parallel
% neuron            complemented neuron structure ready for parallel NEURON
% mindelay          minimum delay and time step of parallel network
%                   (calculated from minimal delay between synapses)
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


mindelay = 10;  % default minimum delay and time step of parallel network
if ~exist('thesetrees','var') || isempty(thesetrees)
    thesetrees = 1:numel(tree);
end
% give each existing cell a unique GID and initialize some basic parameters
% this is necessary that gid_exist works
for t = 1:numel(thesetrees)
    if isfield(tree{thesetrees(t)},'artificial')
        GIDs(t) = struct('gid',t-1,'pp',[],'watch','on','node',[],'cell',t);
    else
        GIDs(t) = struct('gid',t-1,'pp',[],'watch','v','node',1,'cell',t);
    end
end
counter = 1+t; 

if isfield(neuron,'con')
    for c = 1:numel(neuron.con)
        if isfield(neuron.con(c).source,'pp') && ~isempty(neuron.con(c).source.pp)
            % check for all fields and put standard values in there if not
            % existing
            if isfield(neuron{n}.con(c).source,'ppg')  % check for an index to a PP subgroup
                ppg = neuron{n}.con(c).source.ppg;
            else
                ppg = 1:numel(neuron{x}.pp{neuron.con(c).source.cell}.(neuron.con(c).source.pp));  % else take all PP subgroups
            end
            % add the point process which is used as a netcon as a new GID            
            ind = find(arrayfun(@(x) strcmp(x.pp,neuron.con(c).source.pp) & x.cell  == neuron.con(c).source.cell & isequal(x.node,neuron.con(c).source.node) & strcmp(x.watch,neuron.con(c).source.watch) ,GIDs),1,'first');
            if ~isempty(ind)
                neuron.con(c).source.gid = GIDs(ind).gid(1:numel(ppg));
            else
                for p = 1:numel(ppg)
                    neuron.con(c).source.gid(p) = counter-1;
                    GIDs(counter) = neuron.con(c).source;
                    GIDs(counter).gid = counter -1;  % only one gid there..
                    counter = counter +1;
                end
            end
        else
            % check if not pp, check if GID of cell exist, check if it is
            % about the same node, check if the same variable is watched
            ind = arrayfun(@(x) isempty(x.pp) & x.cell  == neuron.con(c).source.cell & isequal(x.node,neuron.con(c).source.node) & strcmp(x.watch,neuron.con(c).source.watch) ,GIDs);
            if any(ind)
                neuron.con(c).source.gid = GIDs(ind).gid;
            else
                neuron.con(c).source.gid = counter-1;
                GIDs(counter) = neuron.con(c).source;
                counter = counter +1;
            end
        end
        % check if there is any netcon delay smaller than mindelay    
        if isfield(neuron.con(c),'delay')
            mindelay = min(mindelay,neuron.con(c).delay);
        end
    end
end
if mindelay <= 0
    error('The minimum delay of NetCons in parallel NEURON has to be greater than zero!')
end
[~,ind] = sort(cat(1,GIDs.cell));
GIDs = GIDs(ind);

end

