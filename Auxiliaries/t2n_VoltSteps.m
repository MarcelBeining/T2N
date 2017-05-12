function t2n_VoltSteps(neuron,tree,params,targetfolder_data,ostruct)
% This function performs one or multiple voltage steps in the cells given
% by "tree" and "neuron" and saves the results in a mat file named
% according to neuron.experiment
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


if nargin < 5 || ~isfield(ostruct,'holding_voltage')
    ostruct.holding_voltage = -80;
end
elecnode = 1;
if ~isfield(ostruct,'dur')
    ostruct.dur = [105 100 105];
end
params.tstop = sum(ostruct.dur);

LJP = params.LJP;  % use LJP only for current clamps now...
if ostruct.coarse
    vstepsModel =  (-130:10:-40) - LJP;
else
    vstepsModel =  (-130:5:-40) - LJP;
end
holding_voltage = ostruct.holding_voltage - LJP;

nneuron = cell(numel(vstepsModel),1);
for s = 1:numel(vstepsModel)
    if s == 1
        nneuron{s} = neuron;
    end
    amp = cat(2,holding_voltage, vstepsModel(s), holding_voltage);
    for t = 1:numel(tree)
        nneuron{s}.pp{t}.SEClamp = struct('node',elecnode,'rs',15,'ostruct.dur', ostruct.dur,'amp', amp);
        nneuron{s}.record{t}.SEClamp = struct('record','i','node',elecnode);
    end
    nneuron{s} = t2n_as(1,nneuron{s});
end
out = t2n(tree,params,nneuron,'-q-d-w');
if any(cellfun(@(x) x.error,out(cellfun(@(x) isfield(x,'error'),out))))
    return
end

ind = vstepsModel == holding_voltage;
for t = 1:numel(tree)
    mholding_current(t) = mean(out{ind}.record{t}.SEClamp.i{1}(find(out{ind}.t>=180,1,'first'):find(out{ind}.t>=200,1,'first')) *1000); % + 6 * mean(outleaksub{s}.record{t}.i{1}(find(outleaksub{s}.t>=180,1,'first'):find(outleaksub{s}.t>=200,1,'first')) *1000);
end
for s = 1:numel(vstepsModel)
    for t = 1:numel(tree)
        steadyStateCurrVec(s,t) =  mean(out{s}.record{t}.SEClamp.i{1}(find(out{s}.t>=180,1,'first'):find(out{s}.t>=200,1,'first')) *1000);% - mholding_current(t); % + 6 * mean(outleaksub{s}.record{t}.i{1}(find(outleaksub{s}.t>=180,1,'first'):find(outleaksub{s}.t>=200,1,'first')) *1000);
        currVec{t,s} =  [out{s}.t';out{s}.record{t}.SEClamp.i{1}' *1000];% - mholding_current(t); % leak subtraction not done by Mongiat + 6 * [0*outleaksub{s}.t';mask.* outleaksub{s}.record{t}.i{1}' *1000];
    end
end
save(fullfile(targetfolder_data,sprintf('Exp_VoltSteps_%s.mat',neuron.experiment)),'mholding_current','neuron','holding_voltage','steadyStateCurrVec','currVec','params','vstepsModel','tree','LJP')
