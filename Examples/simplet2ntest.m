% this script is to show you some applications of t2n. the first two sections initialize the morphologie and parameters, 
% whereas the other sections make different simulations. Run the second
% sections best each time before you run a simulation.

%% load a tree

tree = {sample_tree};   % load sample tree
tree{1}.R(1:2) = 3;     % make the first two nodes a new region
tree{1}.rnames{3} = 'soma';   % name the region "soma"
tree{1}.D(1:2) = 10;        % increase diameter of the soma

figure;plot_tree(tree{1},tree{1}.R)   % plot tree with region code
axis off

%% initialize parameters and neuron structure
params = [];            % clear params
params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';   % add path to NEURON
params.morphfolder = 'morphos/testmorphs';    % where should the transformed trees be saved to?
params.v_init = -80;                    % starting voltage of all cells
params.dt = 0.025;                      % integration time step [ms]
params.tstop = 300;                     % stop simulation after this (simulation) time [ms]
params.prerun = -400;                   % add a pre runtime [ms] to let system settle
params.celsius = 37;
params.nseg = 'd_lambda';
params.accuracy = 1;

neuron = [];                    % clear neuron
for t = 1:numel(tree)
    neuron.mech{t}.all.pas = struct('g',0.0001,'Ra',200,'e',-80,'cm',1);    % add passive channel to all regions and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    neuron.mech{t}.soma.hh = struct('gnabar',0.4,'gl',0);                   % add Hodgkin-Huxley Sodium channel only to soma
end

neuron.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[0.2 0]);  % add a current clamp electrode to first node and define stimulation times [ms] and amplitude [nA]
% neuron.pp{1}.AlphaSyn = struct('node',16,'gmax',0.2,'onset',50);
neuron.record{1}.cell = struct('node',1,'record','v');                  % record voltage "v" from first node
neuron.record{1}.IClamp = struct('node',1,'record','i');

tree = t2n_writetrees(params,tree);                                     % transform tree to NEURON morphology

%% execute neuron and plot
out = t2n(tree,params,neuron,'-w-q');         % execute t2n

figure; 
subplot(2,1,1)
plot(out.t,out.record{1}.cell.v{1})         % plot result (time vs voltage)
ylabel('Voltage [mV]')
xlabel('Time [ms]')
subplot(2,1,2)
plot(out.t,out.record{1}.IClamp.i{1})         % plot result (time vs voltage)

%% do several simulations with different current input
%% t2n CAN DO PARALLEL! :-D

amp = 0:0.04:0.2;       % define a series of amplitudes for the electrode
nneuron = cell(0);      % initialize simulation list
for n = 1:numel(amp)
    nneuron{n} = neuron;        % use previously defined neuron strucutre as basis for each simulation list entry
    nneuron{n}.pp{1}.IClamp.amp = [amp(n),0];   % only change the amplitude of the already inserted current clamp electrode
end

out = t2n(tree,params,nneuron,'-w-q');  % execute t2n

figure; hold all
for n = 1:numel(amp)       % go through simulations/amplitudes
    subplot(numel(amp),1,n)     % make subplot for each simulation
    plot(out{n}.t,out{n}.record{1}.cell.v{1})  % plot time vs voltage of current simulation
    legend(sprintf('%g nA current injection',amp(n)))   % add legend
    ylabel('Voltage [mV]')
end
linkaxes    % make all axes the same size
xlabel('Time [ms]')

%% Map voltage during spike onto tree

t = [50 60];
neuron.record{1}.cell.node = 1:numel(tree{1}.X);  % again use already defined neuron structure but now record from all nodes of the tree

out = t2n(tree,params,neuron,'-w-q');   % execute t2n

% [~, ind] = max(out.record{1}.cell.v{1});   % find the time of the spike maximum at the soma
for tt = 1:numel(t)
    figure;
    plot_tree(tree{1},cellfun(@(x) x(out.t==t(tt)),out.record{1}.cell.v))  % plot the tree with the voltage at that time at each node
    colorbar                % add a colorbar
    set(gca,'CLim',[-50 40])  % make the colors go only from 0 to 50 mV
    axis off
end