function t2n_bAP(neuron,tree,params,targetfolder_data,ostruct)

params.v_init = -85.4;
if ~isfield(ostruct,'cstep')
    ostruct.cstep = 1.3; % nA
end

params.tstop = 1000;
params.dt=0.025;
params.cvode = 1;
nodes = cell(numel(tree),1);
plen = nodes;
eucl = nodes;
plotvals = nodes;  
bAP = nodes;
ipar = nodes;

hstep = t2n_findCurr(params,neuron,tree,params.v_init); %assuming a HP of xxx mV


for t = 1:numel(tree)
    
    plen{t} = Pvec_tree(tree{t});
    ipar{t} = ipar_tree(tree{t});
    ipar{t} = ipar{t}(T_tree(tree{t}),:);  % only paths from termination points
    ipar{t}(ipar{t}==0) = 1;
    
    if ostruct.simple
        nodes{t} = unique(ipar(1,:));
    else
        nodes{t} = unique(ipar{t});
    end
    
    % this part would have been done by t2n anyway, however to avoid
    % loading a lot of redundant values into Matlab, nodes are reduced to
    % the locations were NEURON actually calculates voltage here
    minterf = load(fullfile(params.path,params.morphfolder,sprintf('%s_minterf.mat',tree{t}.NID)));
    minterf = t2n_make_nseg(tree{t},minterf.minterf,params,neuron.mech{t});
    inode = zeros(numel(nodes{t}),1);
    for in = 1:numel(nodes{t})
        inode(in) = find(minterf(:,1) == nodes{t}(in),1,'first');    %find the index of the node in minterf
    end
    [~,ia] = unique(minterf(inode,[2,4]),'rows');
    nodes{t} = sort(nodes{t}(ia));
    if ostruct.reduce  % reduce number of real recorded nodes to every third.
        nodes{t} = nodes{t}(1:3:end);
    end
    
    neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    neuron.pp{t}.IClamp = struct('node',1,'times',[-200 30,32.5],'amp', [hstep(t) hstep(t)+ostruct.cstep hstep(t)]); %n,del,dur,amp
    eucl{t} = eucl_tree(tree{t});
end
[out, ~] = t2n(tree,params,neuron,'-w-q-d');
tim = out.t;
for t = 1:numel(tree)
    plotvals{t} = NaN(numel(tree{t}.X),1);
    for x = 1:numel(nodes{t})
        [mx, ind] = max(out.record{t}.cell.v{nodes{t}(x)});
        basl = mean(out.record{t}.cell.v{nodes{t}(x)}(tim>=0 & tim<=30));
        ind2 = find(out.record{t}.cell.v{nodes{t}(x)} > (mx-basl)/2 + basl,1,'first');
        
        bAP{t}(x,:) = [nodes{t}(x) tim(ind) plen{t}(nodes{t}(x)) eucl{t}(nodes{t}(x)) mx basl tim(ind2)]; % ind nodes, time of max amp, PL at nodes, eucl at nodes, max amplitude, baseline, time of halfmax amp
        plotvals{t}(nodes{t}(x)) = mx;
    end
    % as not all nodes were recorded, interpolate the value for the nodes
    % between the recorded ones (only affects plotting the tree, not the
    % data graphs)
    for x = 1:size(ipar{t},1)
        plotvals{t}(ipar{t}(x,:)) = interp1(plen{t}(intersect(nodes{t},ipar{t}(x,:))),plotvals{t}(intersect(nodes{t},ipar{t}(x,:))),plen{t}(ipar{t}(x,:)),'pchip');
    end
    
end

save(fullfile(targetfolder_data,sprintf('Exp_bAP_%s.mat',neuron.experiment)),'bAP','plotvals','params','nodes','neuron','tree','tim')
