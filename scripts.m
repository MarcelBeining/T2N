
[tree, name, path] = load_tree;
tree = tree{1};
t= 1;

%%
figure
plot_tree(tree{t})

%% Number Branch Points

sum(B_tree(tree{t}))
%%
sum(B_tree(tree{t}) & tree{t}.R == 3)

%% Branch Order
BO_tree(tree{t},'-s')


%% Hull
figure;
hull_tree(tree{t}, 50,[],[],[],'-s')
hold on, plot_tree(tree{t});

%% Steady-State signature
figure;
p = 1021;
sse_tree(tree{t},p,'-s')
hold on
pointer_tree(tree{t},p,[],[],[],'-l')


%%     PRESENTATION

%%
for i = 1:numel(tree)
    tree{t} = sort_tree(tree{t},'-LO');
    tree{t}.R(1:20) = 7;
    tree{t}.rnames = cat(2,tree{t}.rnames,{'soma'});
end
    
%% Export to neuron
del = 300;     %ms
dur = 600;     %ms
amp = 0.015;   %nA

params.prerun = 0;
params.v_init = -76;
params.tstart = 0;
params.tstop = 1000;
params.dt = 0.025;
params.nseg = 'd_lambda';
params.openNeuron = 1;
params.morphfolder = 'morphos/';
params.exchfolder = 't2nexchange_IMPRS';
params.neuronpath = 'C:\nrn73w64\bin64\nrniv.exe';
neuron.pas{1} = [1,80,1/50000,-76.6];

% AH Na+K channels
neuron.act{1}.soma.ichan3 = {'gnabar',0.12,'gkfbar',0.016,'gksbar',0.003,'gkabar',0.012};
neuron.act{1}.GCL.ichan3 = {'gnabar',0.018,'gkfbar',0.004,'gksbar',0.003,'gkabar',0};
neuron.act{1}.SGCL.ichan3 = {'gnabar',0.018,'gkfbar',0.004,'gksbar',0.003,'gkabar',0};
neuron.act{1}.bdend.ichan3 = {'gnabar',0.018,'gkfbar',0.004,'gksbar',0.003,'gkabar',0};
% neuron.act{1}.axon.ichan3 = {'gnabar',0.21,'gkfbar',0.028,'gksbar',0,'gkabar',0.004};
neuron.act{1}.IML.ichan3 = {'gnabar',0.13,'gkfbar',0.004,'gksbar',0.003,'gkabar',0};
neuron.act{1}.MML.ichan3 = {'gnabar',0.08,'gkfbar',0.001,'gksbar',0.003,'gkabar',0};
neuron.act{1}.OML.ichan3 = {'gnabar',0,'gkfbar',0.001,'gksbar',0.004,'gkabar',0};
neuron.act{1}.OMLout.ichan3 = {'gnabar',0,'gkfbar',0.001,'gksbar',0.004,'gkabar',0};

%AH Calcium channels
neuron.act{1}.soma.Ca = {'gtcabar', 0.00015, 'gncabar', 0.002,'glcabar', 0.01};
neuron.act{1}.SGCL.Ca = {'gtcabar', 0.0003, 'gncabar', 0.003,'glcabar', 0.015};
neuron.act{1}.GCL.Ca = {'gtcabar', 0.0003, 'gncabar', 0.003,'glcabar', 0.015};
neuron.act{1}.bdend.Ca = {'gtcabar', 0.0003, 'gncabar', 0.003,'glcabar', 0.015};
neuron.act{1}.IML.Ca = {'gtcabar', 0.001, 'gncabar', 0.001,'glcabar', 0.015};
neuron.act{1}.MML.Ca = {'gtcabar', 0.002, 'gncabar', 0.001,'glcabar', 0.001};
neuron.act{1}.OML.Ca = {'gtcabar', 0.002, 'gncabar', 0.001,'glcabar', 0};

neuron.act{1}.soma.CadepK = {'gbkbar', 0.0003,'gskbar', 0.0005};
neuron.act{1}.SGCL.CadepK = {'gbkbar', 0.0003,'gskbar', 0.0002};
neuron.act{1}.GCL.CadepK = {'gbkbar', 0.0003,'gskbar', 0.0002};
neuron.act{1}.bdend.CadepK = {'gbkbar', 0.0003,'gskbar', 0.0005};
neuron.act{1}.IML.CadepK = {'gbkbar', 0.0005,'gskbar', 0.0001};
neuron.act{1}.MML.CadepK = {'gbkbar', 0.0012,'gskbar', 0};
neuron.act{1}.OML.CadepK = {'gbkbar', 0.0012,'gskbar', 0};

neuron.act{1}.na_ion = {'ena',45};
neuron.act{1}.k_ion = {'ek',-85};

neuron.stim{1} = [1, del, dur, amp];


[out, minterf] = t2n(tree{t},params,neuron,path,'-d');
