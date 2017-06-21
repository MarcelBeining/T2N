% this script is to show you some applications of t2n. 
% the first two sections initialize the morphology and standard parameters, 
% whereas the other sections define the different simulations. 

%% load a tree
% First we need a morphology (tree), on which we can add the channels etc.
% We are using here the sample_tree function from the TREES toolbox which
% loads a small tree. You can also load a different tree by using
% tree = load_tree

tree = {sample_tree};   % load sample tree
% this tree has no soma regino yet, so just make the first two nodes
% "somatic"
tree{1}.R(1:2) = 3;     % make the first two nodes a new region
tree{1}.rnames{3} = 'soma';   % name the region "soma"
tree{1}.D(1:2) = 10;        % increase diameter of the somatic nodes

figure;plot_tree(tree{1},tree{1}.R),colorbar   % plot tree with region code
axis off


%% initialize parameters and neuron structure
% now we need to initialize the t2n parameter structure which includes
% general NEURON settings. Most of these settings are set to a default
% value if not explicitly set, but you might want to control most of them

params = [];                                         % clear params structure
params.neuronpath = 'C:/nrn73w64/bin64/nrniv.exe';   % add path to NEURON exe (only necessary for Windows)
params.morphfolder = 'morphos/testmorphs';           % relative path of where the transformed trees are saved to
params.v_init = -80;                                 % starting membrane potential [mV] of all cells
params.dt = 0.025;                                   % integration time step [ms]
params.tstop = 300;                                  % stop simulation after this (simulation) time [ms]
params.prerun = -400;                                % add a pre runtime [ms] to let system settle
params.celsius = 10;                                 % temperature [celsius]
params.nseg = 'dlambda';                             % the dlambda rule is used to set the amount of segments per section. Alternatively, a fixed number can be entered or a string 'EachX' with X being the interval between segments [micron]
params.accuracy = 0;                                 % optional argument if number of segments should be increased in the soma and axon

% now we set the t2n neuron structure which defines all mechanisms, point
% processes, synapses, connections recordings, protocols etc. that should
% be used in the simulation. The mechanisms and point processes can be any
% that is defined in the standard NEURON environment or that has been 
% written as .mod file and saved to the 'lib_mech' folder in the model
% folder from which it can be compiled by t2n. 'neuron' is a Matlab structure 
% if a single simulation should be run 

neuron = [];                                                                % clear neuron structure
for t = 1:numel(tree)                                                       % loop though all morphologies (is one when using sample_tree)
    neuron.mech{t}.all.pas = struct('g',0.0003,'Ra',100,'e',-80,'cm',1);    % add passive channel to all regions and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    % alternatively to the compact code above, one can define each
    % parameter line by line as with NEURON by writing:
    % neuron.mech{t}.all.pas.g = 0.0003;
    % neuron.mech{t}.all.pas.Ra = 100;
    % neuron.mech{t}.all.pas.e = -80;
    % neuron.mech{t}.all.pas.cm = 1;
    neuron.mech{t}.all.k_ion.ek = -90;
    neuron.mech{t}.all.na_ion.ena = 50;
    neuron.mech{t}.soma.hh = struct('gnabar',0.25,'gkbar',0.036,'gl',0);     % add Hodgkin-Huxley Sodium channel only to soma
    neuron.record{t}.cell = struct('node',1,'record','v');                   % record voltage "v" from first node (i.e. the soma)
end

tree = tree(1);                                                              % for simplicity, only one tree is looked at (if several exist at all)

tree = t2n_writeTrees(tree,params,fullfile(pwd,'test.mtr'));                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology


%% Tutorial 1 - simulation protocol: somatic current injection
% Now we want to do a simple somatic current injection, execute neuron and 
% plot the result. In order that we do not have to rerun the upper sections
% each time, we will just copy the 'neuron' structure into a new variable
% 'nneuron' and modify only this one in each protocol.
nneuron = neuron;                                                           % copy standard neuron structure

% now we define a 100 ms current injection of 0.6 nA that is started at 
% time point 50 ms
nneuron.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[0.6 0]);     % add a current clamp electrode to the first node and define stimulation times [ms] and amplitudes [nA]
nneuron.record{1}.IClamp = struct('node',1,'record','i');                   % record the current of the IClamp just for visualization

out = t2n(tree,params,nneuron,'-w-q');                                      % execute t2n and receive output

% After execution of t2n we can plot the results. T2N returns all
% recordings that had previously been defined. We can access them in
% out.record{1} where {1} accesses the recordings of the first
% tree/morphology. Subsequently all recordings of distributed mechanisms or
% other cell-specific parameters (e.g. voltage) can be found in the field
% 'cell', whereas recorded values from point processes can be found under
% the name of that process (e.g. 'IClamp'). Following that comes the name
% of the recorded variable and the index to the node at which the variable
% has been recorded or the corresponding point process had been placed
% (e.g. v{1} for the voltage at the first node). The time vector can be
% found in out.t and fits to all recordings, except if the options
% params.local_dt has been used which allows NEURON to use different time
% steps for different trees. In that case out.t is a cell array of vectors.

figure; 
subplot(2,1,1)                                                              % make a subplot in the figure
plot(out.t,out.record{1}.cell.v{1})                                         % plot recorded somatic voltage (time vs voltage)
ylim([-90,50])
xlim([0,params.tstop])
ylabel('Membrane potential [mV]')
xlabel('Time [ms]')
subplot(2,1,2)                                                              % make another subplot in the figure
plot(out.t,out.record{1}.IClamp.i{1})                                       % plot electrode current (time vs current)
ylim([0,1])
xlim([0,params.tstop])
ylabel('Injected current [nA]')


%% Tutorial 2 - simulation protocol: Do several simulations with different injected current amplitudes
% Now we learn something about the real strength of t2n: It can execute 
% simulations in parallel! All we have have to do is to create a cell array
% of neuron structures (each defining a different protocol) and hand them 
% to t2n. Here we use this to simulate different current injections, but in
% principle this can be used for any different protocols (e.g. different
% channels, different synapses, different stimulation patterns etc).

amp = 0:0.12:0.6;                  % define a series of amplitudes for the electrode (here: 6 steps from 0 to 0.6 nA)
nneuron = cell(numel(amp),1);      % initialize (clear) simulation list

% principally, we could copy the standard neuron structure into the neuron cell 
% array for each protocol and then add the IClamp with different amplitudes
% which would like this: 
% for n = 1:numel(amp)
%     nneuron{n} = neuron;                                                            % use neuron structure defined above as basis
%     nneuron{n}.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[amp(n) 0]);   % only define what is changed in for the next simulations
% end
% However, in this case T2N would translate the complete descriptions of
% all mechanisms etc into hoc code for each simulation. This might be no
% big deal when defining only 6 simulations, but could be very time and
% (to lesser extent) disk space consuming when defining hundreds of 
% simulations. Hence, it is recommended to use the t2n_as function which
% tells t2n to simply reuse the hoc file of a specific simulation instance
% for all definitions that have not been defined in the other simulations.
% Here is how it is used:

nneuron{1} = neuron;                                                               % use neuron structure defined above as basis for the first simulation instance
for n = 1:numel(amp)
    nneuron{n}.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[amp(n) 0]);   % now define the IClamp in all simulations (remember that all other parameters were only defined in simulation 1) and...
end
nneuron = t2n_as(1,nneuron);                                                        % ... use the previously defined neuron structure as basis for each other simulation

% now T2N only writes hoc code for initialization of the IClamp point 
% process in each simulation and uses the hoc code about morphology, 
% recordings and mechanisms from the first simulation

out = t2n(tree,params,nneuron,'-w-q');                                              % execute t2n

% Now we can plot the recorded somatic voltages for each current step
figure; hold all
for n = 1:numel(amp)                                    % go through simulations/amplitudes
    subplot(numel(amp),1,n)                             % make subplot for each simulation
    plot(out{n}.t,out{n}.record{1}.cell.v{1})           % plot time vs voltage of current simulation
    title(sprintf('%g nA current injection',amp(n)))    % add legend
    ylabel('Membrane pot. [mV]')
end
linkaxes                                                % make all axes the same size
xlabel('Time [ms]')

%% Tutorial 3 - Map the membrane potential during a spike onto the tree
% As t2n features the trees toolbox, morphological visualizations are easy.
% We illustrate this by mapping the membrane potential at specific times 
% during a current injection on our morphology.

nneuron = neuron;                                                         % copy standard neuron structure
nneuron.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[0.6 0]);   % define the IClamp
times = [52 53 54 55];                                                    % define some time points at which the voltage should be mapped
nneuron.record{1}.cell.node = 1:numel(tree{1}.X);                         % use already defined recording but modify it to record from all nodes of the tree

out = t2n(tree,params,nneuron,'-w-q');                                    % execute t2n

% The plot_tree function of the TREES toolbox allows handing over a vector
% with the same size as the number of nodes which is then mapped on the
% tree. This means we need to get only one value from each recorded voltage
% vector. You can do this by using the Matlab cellfun function, or simply
% use our t2n_get function which can be used to extract values at a given
% time (or time span) or to easily apply a function (such as 'mean') on
% each recorded vector (see documentation of t2n_get).
figure;
for tt = 1:numel(times)                                                   % loop through all defined time points
    subplot(ceil(sqrt(numel(times))),ceil(sqrt(numel(times))),tt)         % make quadratically arranged subplots
    vec = t2n_get(out,'v',times(tt));                                     % extract the voltages at that time point
    plot_tree(tree{1},vec{1})                                             % plot the tree with the voltage at that time at each node
    colorbar                                                              % add a colorbar
    set(gca,'CLim',[-50 40])                                              % make the colors go only from -50 to 40 mV
    axis off
    title(sprintf('Tree at time %g ms',times(tt)))
end
%% Tutorial 4 - parameter scan
% Now we want to analyze the impact of a parameter on spiking, e.g. for
% optimizing the model. We simply use a for loop to create neuron
% instances, each being the same except for one parameter that is changed.

fac = 0:0.5:2;                                          % define the multiplication factor range of the parameter scan
nneuron = cell(numel(fac),1);                           %initialize neuron cell array (generally good strategy if you are using this variable in different sections)
for f = 1:numel(fac)                                    % loop through number of different factors
    nneuron{f} = neuron;                                % copy standard neuron definition
    nneuron{f}.mech{1}.soma.hh.gkbar = fac(f) * 0.036;     % change the HH potassium channel density according to fac
end
nneuron{1}.pp{1}.IClamp = struct('node',1,'times',[50 150],'amp',[0.6 0]);    % again define the IClamp only for first instance...
nneuron = t2n_as(1,nneuron);                                                  % ... and let t2n use it for the rest

out = t2n(tree,params,nneuron,'-w-q');   % execute t2n

% plot the result
figure; hold all
for f = 1:numel(fac)
    plot(out{f}.t,out{f}.record{1}.cell.v{1})  % plot time vs voltage of each simulation
end
legend(arrayfun(@(x) sprintf('hh K-channel factor of %g',x),fac,'uni',0))   % add legend
xlabel('Time [ms]')
ylabel('Voltage [mV]')

%% Tutorial 5 - Synaptic stimulation, simple Alpha synapse
% Now we come closer to the point of building networks. First we do a
% synaptic stimulation of a chosen point within the morphology using the
% NEURON AlphaSynapse point process which adds synaptic current at a 
% specific time "onset" at the location it is placed. For this example, we
% chose the most far away point in the tree as the synapse location, which 
% can be found by using the TREES toolbox function Pvec_tree. Other
% examples to extract a specific synapse location that can be tested here
% is: searching for termination points with "find(T_tree(tree{1})"; 
% getting nodes at a specific distance to the soma with "find(plen == x)" 
% (or add some tolerance y: "find(plen-x < y)"); 
% getting nodes of a specific region name REG with "find(tree{1}.R ==
% find(strcmp(tree{1}.rnames,REG))" or any combinaton of them.

nneuron = neuron;                                                         % copy standard neuron structure
plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy    
nneuron.pp{1}.AlphaSynapse = struct('node',synIDs,'gmax',0.01,'onset',20);% add an AlphaSynapse at this location
nneuron.record{1}.cell.node = cat(1,1,synIDs);                            % record somatic voltage and voltage at synapse
nneuron.record{1}.AlphaSynapse = struct('node',synIDs,'record','i');      % record synaptic current

out = t2n(tree,params,nneuron,'-w-q');              % execute t2n

% plot the result (Vmem at soma and synapse and synaptic current)
figure;
subplot(2,1,1)
hold all
plot(out.t,out.record{1}.cell.v{1})                 % plot time vs voltage at soma
plot(out.t,out.record{1}.cell.v{synIDs})            % plot time vs voltage at dendrite end
legend('Soma','Synapse')
ylabel('Membrane potential [mV]')
ylim([-85,-60])
xlim([0,50])
subplot(2,1,2)
plot(out.t,out.record{1}.AlphaSynapse.i{synIDs})    % plot time vs synaptic current
xlim([0,50])
ylabel('Synaptic current [nA]')
xlabel('Time [ms]')

%% Tutorial 6 - Synaptic stimulation, Exp2Syn synapse and a NetStim
% A more customizable synapse in NEURON is the Exp2Syn point process as it
% is activated by another cell or point process. Hence, if we use such a
% synapse, we have to make use of NEURON NetCon objects, which can be
% defined in the .con field of the neuron structure. Additionally, we have
% to introduce some cell that stimulates the synapse at given times, for
% which we will introduce a NetStim as an artificial cell. This cell, as
% all real and artificial cells, has to be defined in the tree structure
% and is quite simple to define:

tree{2} = struct('artificial','NetStim','params',struct('start',10,'interval',15,'number',10)); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
tree = t2n_writeTrees(tree,params,fullfile(pwd,'test.mtr'));                                    % tree morphologies are rewritten because this NetStim might be not written yet

nneuron = neuron;                                                         % copy standard neuron structure
plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy    
nneuron.pp{1}.Exp2Syn = struct('node',synIDs,'tau1',0.2,'tau2',2.5,'e',0);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
nneuron.record{1}.cell.node = cat(1,1,synIDs);                            % record somatic voltage and voltage at synapse
nneuron.record{1}.Exp2Syn = struct('node',synIDs,'record','i');           % record synaptic current
nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs),'delay',0,'threshold',0.5,'weight',0.1);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)

out = t2n(tree,params,nneuron,'-w-q');                                    % execute t2n

% plot the result (Vmem at soma and synapse and synaptic current)
figure;
subplot(2,1,1)
hold all
plot(out.t,out.record{1}.cell.v{1})       % plot time vs voltage at soma
plot(out.t,out.record{1}.cell.v{synIDs})  % plot time vs voltage at dendrite end
legend('Soma','Synapse')
ylabel('Membrane potential [mV]')
xlabel('Time [ms]')
% ylim([-85,-60])
% xlim([0,100])
subplot(2,1,2)
plot(out.t,out.record{1}.Exp2Syn.i{synIDs})  % plot time vs synaptic current


%% Tutorial 7 - Synaptic stimulation, Exp2Syn synapse and a NetStim
% A more customizable synapse in NEURON is the Exp2Syn point process as it
% is activated by another cell or point process. Hence, if we use such a
% synapse, we have to make use of NEURON NetCon objects, which can be
% defined in the .con field of the neuron structure. Additionally, we have
% to introduce some cell that stimulates the synapse at given times, for
% which we will introduce a NetStim as an artificial cell. This cell, as
% all real and artificial cells, has to be defined in the tree structure
% and is quite simple to define:

tree{2} = struct('artificial','NetStim','params',struct('start',10,'interval',15,'number',10)); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
tree = t2n_writeTrees(tree,params,fullfile(pwd,'test.mtr'));                                    % tree morphologies are rewritten because this NetStim might be not written yet

nneuron = neuron;                                                         % copy standard neuron structure
plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy    
nneuron.pp{1}.Exp2Syn = struct('node',synIDs,'tau1',0.2,'tau2',2.5,'e',0);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
nneuron.record{1}.cell.node = cat(1,1,synIDs);                            % record somatic voltage and voltage at synapse
nneuron.record{1}.Exp2Syn = struct('node',synIDs,'record','i');           % record synaptic current
nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',synIDs),'delay',0,'threshold',0.5,'weight',0.1);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in synIDs), and add threshold/weight and delay of the connection (NetStim parameters)

out = t2n(tree,params,nneuron,'-w-q');                                    % execute t2n

% plot the result (Vmem at soma and synapse and synaptic current)
figure;
subplot(2,1,1)
hold all
plot(out.t,out.record{1}.cell.v{1})       % plot time vs voltage at soma
plot(out.t,out.record{1}.cell.v{synIDs})  % plot time vs voltage at dendrite end
legend('Soma','Synapse')
ylabel('Membrane potential [mV]')
xlabel('Time [ms]')
% ylim([-85,-60])
% xlim([0,100])
subplot(2,1,2)
plot(out.t,out.record{1}.Exp2Syn.i{synIDs})  % plot time vs synaptic current
