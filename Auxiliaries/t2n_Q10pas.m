function neuron = t2n_Q10pas(neuron,celsius)
% changes passive parameters g and Ra to adjust it to a given temperature
% (assuming that g and Ra were defined for 24C before). Do not use it
% multiple times

scaleg = 1.98^((celsius-24)/10);
scaleRa = 0.8^((celsius-24)/10);

for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    for f1 = 1:numel(fields)
        if isfield(neuron.mech{t}.(fields{f1}),'pas')
            neuron.mech{t}.(fields{f1}).pas.g = neuron.mech{t}.(fields{f1}).pas.g * scaleg;
            neuron.mech{t}.(fields{f1}).pas.Ra = neuron.mech{t}.(fields{f1}).pas.Ra * scaleRa;
        end
    end
end

if isfield(neuron,'experiment') && ~isempty(neuron.experiment)
    neuron.experiment = strcat(neuron.experiment,'_Q10pas');
else
    neuron.experiment = 'Q10pas';
end
