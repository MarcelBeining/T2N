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
for t = 1:numel(thesetrees)
    GIDs(t) = struct('gid',t-1,'pp',[],'watch','v','node',1,'cell',thesetrees(t));
end
counter = 1+t; 

if isfield(neuron,'con')
    for c = 1:numel(neuron.con)
        if isfield(neuron.con(c).source,'pp') 
            if ~isfield(neuron.con(c).source,'watch')
                neuron.con(c).source.watch = 'on';
            end
            if ~isfield(neuron.con(c).source,'node')
                neuron.con(c).source.node = 1;
            end
            if isfield(neuron{n}.con(c).source,'ppg')  % check for an index to a PP subgroup
                ppg = neuron{n}.con(c).source.ppg;
            else
                ppg = 1:numel(neuron{x}.pp{neuron.con(c).source.cell}.(neuron.con(c).source.pp));  % else take all PP subgroups
            end
                        
            ind = find(strcmp({GIDs.pp},neuron.con(c).source.pp) & cat(1,GIDs.cell) == neuron.con(c).source.cell & cat(1,GIDs.node) == neuron.con(c).source.node & strcmp({GIDs.watch},neuron.con(c).source.watch),1,'first');
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
            if ~isfield(neuron.con(c).source,'watch')
                if isfield(tree(thesetrees==neuron.con(c).source.cell),'artificial')
                    neuron.con(c).source.watch = 'on';
                else
                    neuron.con(c).source.watch = 'v';
                end
            end
            if ~isfield(neuron.con(c).source,'node')
                neuron.con(c).source.node = 1;
            end
            ind = cellfun(@isempty,{GIDs.pp}) & cat(1,GIDs.cell) == neuron.con(c).source.cell & cat(1,GIDs.node) == neuron.con(c).source.node & strcmp({GIDs.watch},neuron.con(c).source.watch);
            if any(ind)
                neuron.con(c).source.gid = GIDs(ind).gid;
            else
                neuron.con(c).source.gid = counter-1;
                GIDs(counter) = neuron.con(c).source;
                counter = counter +1;
            end
        end
            
        if isfield(neuron.con(c),'delay')
            mindelay = min(mindelay,neuron.con(c).delay);
        else
            mindelay = min(mindelay,1); % default delay
        end
    end
end

[~,ind] = sort(cat(1,GIDs.cell));
GIDs = GIDs(ind);

end

