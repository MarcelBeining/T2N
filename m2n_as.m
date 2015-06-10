function thisneuron = m2n_as(as,neuron)

fields = {'tree','mech','stim','pp','con','record','play','APCount'};

if nargin < 1 || isempty(neuron)
    thesefields = fields;
else
    thesefields = setdiff(fields,fieldnames(neuron));
    thisneuron = neuron;
end

for f = 1:numel(thesefields)
    if strcmp(thesefields{f},'tree')
        thisneuron.(thesefields{f}) = sprintf('sim%d',as);
    else
        thisneuron.(thesefields{f}) = as;
    end
end