function thisneuron = t2n_as(as,neuron)

fields = {'tree','mech','pp','con','record','play','APCount'};

if nargin < 2 || isempty(neuron)
    thesefields = fields;
else
    thesefields = setdiff(fields,fieldnames(neuron));
    thisneuron = neuron;
end
if nargin < 1
    as = 1;
end

for f = 1:numel(thesefields)
    if strcmp(thesefields{f},'tree')
        thisneuron.(thesefields{f}) = sprintf('sim%d',as);
    else
        thisneuron.(thesefields{f}) = as;
    end
end