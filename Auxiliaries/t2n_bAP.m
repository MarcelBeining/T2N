function t2n_bAP(neuron,tree,params,targetfolder_data,ostruct)

params.v_init = -85.4;
params.skiprun = 0;
if ~(ostruct.vmodel >= 0)  %AH99 model, cstep has to be bigger otherwise not all cells fire
    cstep = 1700*0.001; %nA
else
    cstep = 1300*0.001; %nA
end

params.tstop = 1000;
params.dt=0.025;
params.cvode = 1;
nodes = cell(3,1);
plen = nodes;
eucl = nodes;
plotvals = nodes;  
plotcaivals = nodes;
bAP = nodes;
CaNodes = nodes;
ipar = nodes;
nodes = nodes;

if ~(ostruct.vmodel>=0) 
    cai = 'caim_Caold';  % AH99 model, does not have an own calcium buffer mechanism, thus variable has a different name
    neuron = manipulate_Ra(neuron,1,'axon'); % necessary because otherwise the axon spikes permanently after the buzz and Ca decay measurement is not possible
    params.nseg = 1;  % necessary because dlambda calculation takes to long with high Ra
    params.accuracy = 1; % improve AIS segment number for more accurate simulation
else
    cai = 'cai';
end

hstep = find_curr(params,neuron,tree,params.v_init); %assuming a HP of xxx mV


for t = 1:numel(tree)
    
    plen{t} = Pvec_tree(tree{t});
    CaNodes{t} = {find(abs(plen{t}-150)<0.5 & tree{t}.R == find(strcmp(tree{t}.rnames,'axon'))),find(plen{t} == 0,1,'first'),find(abs(plen{t}-70)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))),find(abs(plen{t}-150)<0.5 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')))};
    if any(cellfun(@isempty,CaNodes{t}))
        fprintf('thresh distance of 0.5 for tree %d did not suffice. switched to thresh of 1\n',t)
        CaNodes{t} = {find(abs(plen{t}-150)<1 & tree{t}.R == find(strcmp(tree{t}.rnames,'axon'))),find(plen{t} == 0,1,'first'),find(abs(plen{t}-70)<1 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon'))),find(abs(plen{t}-150)<1 & tree{t}.R ~= find(strcmp(tree{t}.rnames,'axon')))};
    end
    sprintf('Proximal tree %d: %s',t,strcat(tree{t}.rnames{unique(tree{t}.R(CaNodes{t}{3}))}))
    sprintf('Distal tree %d: %s',t,strcat(tree{t}.rnames{unique(tree{t}.R(CaNodes{t}{4}))}))
    
    ipar{t} = ipar_tree(tree{t});
    ipar{t} = ipar{t}(T_tree(tree{t}),:);  % only paths from termination points
    ipar{t}(ipar{t}==0) = 1;
    
    if ostruct.simple
        uipar = ipar(1,:);
        nodes{t} = uipar(1:find(uipar==1,1,'first'));
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
    nodes{t} = nodes{t}(ia);
    if ostruct.reduce  % reduce number of real recorded nodes to every third.
        nodes{t} = nodes{t}(1:3:end);
    end
    
    neuron.record{t}.cell = struct('node',unique(cat(1,CaNodes{t}{:},nodes{t})),'record',{cai,'v'});
    neuron.pp{t}.IClamp = struct('node',1,'times',[-200 30,32.5],'amp', [hstep(t) hstep(t)+cstep hstep(t)]); %n,del,dur,amp
    eucl{t} = eucl_tree(tree{t});
end
% params.skiprun = 1
[out, ~] = t2n(tree,params,neuron,'-w-q-d');

tw = NaN(numel(tree),4,max(cellfun(@(x) numel(x{4}),CaNodes)));
maxcai = tw;


mCai = cell(numel(tree),numel(CaNodes{1}));
stdCai = mCai;
mV = mCai;
stdV = mCai;

for t = 1:numel(tree)
    plotvals{t} = NaN(numel(tree{t}.X),1);%ones(numel(tree{t}.X),1) * (-80);
    plotcaivals{t} = NaN(numel(tree{t}.X),1);%ones(numel(tree{t}.X),1) * (-80);
    for x = 1:numel(nodes{t})
        [mx, ind] = max(out.record{t}.cell.v{nodes{t}(x)});
        basl = mean(out.record{t}.cell.v{nodes{t}(x)}(out.t>=0 & out.t<=30));
        ind2 = find(out.record{t}.cell.v{nodes{t}(x)} > (mx-basl)/2 + basl,1,'first');
        
        bAP{t}(x,:) = [nodes{t}(x) out.t(ind) plen{t}(nodes{t}(x)) eucl{t}(nodes{t}(x)) mx basl out.t(ind2)]; % ind nodes, time of max amp, PL at nodes, eucl at nodes, max amplitude, baseline, time of halfmax amp
        plotvals{t}(nodes{t}(x)) = mx;
        mx = max(out.record{t}.cell.(cai){nodes{t}(x)});
        plotcaivals{t}(nodes{t}(x)) = mx;
    end
    % as not all nodes were recorded, interpolate the value for the nodes
    % between the recorded ones (only affects plotting the tree, not the
    % data graphs)
    for x = 1:size(ipar{t},1)
        plotvals{t}(ipar{t}(x,:)) = interp1(plen{t}(intersect(nodes{t},ipar{t}(x,:))),plotvals{t}(intersect(nodes{t},ipar{t}(x,:))),plen{t}(ipar{t}(x,:)),'pchip');
        plotcaivals{t}(ipar{t}(x,:)) = interp1(plen{t}(intersect(nodes{t},ipar{t}(x,:))),plotcaivals{t}(intersect(nodes{t},ipar{t}(x,:))),plen{t}(ipar{t}(x,:)),'pchip');
    end
    for f =1:numel(CaNodes{t})
        for ff = 1:numel(CaNodes{t}{f})
            caivec = out.record{t}.cell.(cai){CaNodes{t}{f}(ff)};
            if ~isempty(caivec)  % cai was existent at that node (might not be the case for AH99
%                 [maxcai(t,f,ff),ind(1)] = findpeaks(caivec,'npeaks',1);
%                 % findpeaks does strange things sometimes... use max
                [maxcai(t,f,ff),ind(1)] = max(caivec);
                ind(2) = find(caivec(ind(1):end) <= caivec(end)+ 1e-7,1,'first')+ind(1)-1;
                if diff(ind) > 2
                    est = fit(out.t(ind(1):ind(2))-out.t(ind(1)),caivec(ind(1):ind(2))-caivec(ind(2)),'exp2');
                    %                 est.a/(est.a+est.c)/(-est.b) +  est.c/(est.a+est.c) / (-est.d)
                    if est.a < 0 || est.b < 0
                        est = fit(out.t(ind(1):ind(2))-out.t(ind(1)),caivec(ind(1):ind(2))-caivec(ind(2)),'exp1');
                        tw(t,f,ff) = -1/est.b;
                    else
                        tw(t,f,ff) = (-est.a/est.b+-est.c/est.d)/(est.a+est.c); % amplitude weighted tau as in stocca 2008
                    end
                else
                    tw(t,f,ff) = NaN;
                end
            else
                maxcai(t,f,ff) = NaN;
                tw(t,f,ff) = NaN;
            end
        end
        mCai{t,f} = mean(cat(2,out.record{t}.cell.(cai){CaNodes{t}{f}})*1e6,2)';
        stdCai{t,f} = std (cat(2,out.record{t}.cell.(cai){CaNodes{t}{f}})*1e6,[],2)';
        mV{t,f} = mean(cat(2,out.record{t}.cell.v{CaNodes{t}{f}}),2)';
        stdV{t,f} = std (cat(2,out.record{t}.cell.v{CaNodes{t}{f}}),[],2)';
    end
end
tim = out.t;

save(fullfile(targetfolder_data,sprintf('Exp_bAP_%s.mat',neuron.experiment)),'bAP','plotvals','plotcaivals','params','nodes','neuron','tree','CaNodes','mCai','stdCai','mV','stdV','tim','tw','maxcai')
