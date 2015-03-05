% testm2n
% oldtree = tree;
tree=cell(1);
tree{1}.artificial = 'NetStim';
tree{1}.params.start = 100;
tree{1}.params.interval = 50;
tree{1}.params.number = 100;
tree{2}.artificial = 'NetStim';
tree{2}.params.start = -1000;
tree{2}.params.interval = 5;
tree{2}.params.number = 3;
path = 'D:\Dropbox\PhD Deller\Modeling\P5 GC Model\MarcelCode\yGCs';

tree{3} = oldtree;

neuron.pas{3} = [1, NaN , NaN , -75];
neuron.syn{3} = {'ExpSyn',1,{}};
neuron.syn{3}{1,3}.tau = 0.1;
neuron.con = {'cell',1,'v','cell',2,1,0,5;'cell',2,'v','syn','3.1',1,10,5};
neuron.record{3} = {1,'v'};
neuron.APCount{2} = [1,0];
neuron.APCount{3} = [1,-10];
params = [];
recnode = 2;
params.prerun = 100;
% params.tname = 'GaCa';
params.v_init = -75.3;
params.tstart = 0;
params.tstop = 500;
params.dt = 0.1;
params.nseg = 'd_lambda';
params.openNeuron = 1;
params.exchfolder = 'test';
params.changed = struct('morph',1,'stim',1,'basic',1,'lib',1,'rec',1,'play',1,'connect',1,'act',1,'pas',1,'syn',1,'con',1);


[out, minterf] = m2n(tree,params,neuron,path,'-d');

plot(out.t,out.record{3}.v{1})