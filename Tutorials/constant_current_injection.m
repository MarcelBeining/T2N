% This an example code for the following scenario:
% equip the TREES sample_tree with the HH mechanism in nodes that have a branch order < 2.
% further equip the tree with synapses in all dendritic terminals and inject a current there from
% t = 1sec-2sec. Measure voltage everywhere. 

params = [];
neuron = [];
params.tstop = 3000;
% params.dt = 1;
params.cvode = 1;
params.celsius = 24;
params.nseg = 'd_lambda';
params.exchfolder = 't2nexchange';  % folder were simulation is written
params.morphfolder = 'morphos/NEURON_hocs';  % folder for hoc morphology files

g_pas = 0.0001;  % conductance of passive channel
e_pas = -70;  % e leak
gnabar_hh = 0.0001;  % conductance of hh na (optional)
gkbar_hh = 0.0001;   % conductance of hh k (optional)
amp = 0.1; % amplitude of current injection (nA)

[tree,treename] = load_tree('sample.mtr','none');
tree = resample_tree(tree,1);

if isstruct(tree)
    tree = {tree};
end

%% add mechanisms, point processes etc
for t = 1:numel(tree)
    
   BO = BO_tree(tree{t}) ;
   tree{t}.R(BO < 2) = max(tree{t}.R) + 1;   % give region for hh insertation an own region name (most easy for mechanism insertion)
   tree{t}.rnames = cat(2,tree{t}.rnames,'BOst2');
   
   termind = find(T_tree(tree{t}));  % get terminals
   
   neuron.mech{t}.all.pas = struct('g',g_pas,'e',e_pas,'cm',1);
   neuron.mech{t}.BOst2.hh = struct(); % struct('gnabar',gnabar_hh,'gkbar',gkbar_hh); % wenn du hh nicht mit standardparameter willst ersetze = struct(); mit dem auskommentierten
%    neuron.pp{t}.AlpaSynapse = struct('node',termind);  % this would have been the code for inserting real synapses but for a constant current injection, Iclamp is easier.
   neuron.pp{t}.IClamp = struct('node',termind,'del',1000,'dur',1000,'amp',amp);  % add current injection to terminal nodes
   
   neuron.record{t}.cell = struct('node',1:numel(tree{t}.X),'record','v');
end


%% rewrite tree hocs if necessary ( if morphology changed)
tree = t2n_writetrees(params,tree,fullfile(params.morphfolder,treename));


%% run simulation and plot stuff
out = t2n(tree,params,neuron);

[~,spikeind] = findpeaks(out.record{1}.cell.v{1},'MinPeakHeight',-10);

figure;plot(out.t,out.record{1}.cell.v{1})
xlabel('Time [ms]')
ylabel('Membrane potential [mV]')

figure;plot_tree(tree{1},cellfun(@(x) x(spikeind),out.record{1}.cell.v))
title('Membrane voltage during spike')
colorbar