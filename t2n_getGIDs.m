function  [GIDs,neuron,mindelay] = t2n_getGIDs(neuron,trees,thesetrees)


mindelay = 10;  % default minimum delay and time step of parallel network

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
                if isfield(trees(thesetrees==neuron.con(c).source.cell),'artificial')
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

