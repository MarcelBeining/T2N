function newneuron = t2n_as(x,neuron)
% this function adds all missing fields to a neuron structure by
% telling T2N to use the same fields as those of the xth neuron structure

fields = {'tree','mech','pp','con','record','play','APCount'};

if nargin < 2 || isempty(neuron)
    thesefields = fields;
else
    thesefields = setdiff(fields,fieldnames(neuron));
    newneuron = neuron;
end
if nargin < 1
    x = 1;
end

for f = 1:numel(thesefields)
    if strcmp(thesefields{f},'tree')
        newneuron.(thesefields{f}) = sprintf('sim%d',x);
    else
        newneuron.(thesefields{f}) = x;
    end
end