function strct = changemech(strct,change,mode)
% this function changes parameters in the T2N neuron mech structure according to input change.
% INPUT
% strct: neuron structure containing ion channel densities 
% change: structure with field names according to mechanism parameter names that should be changed and values either describing the factor by which it should be changed or an absolute value
% mode: how parameters should be changed: 1 or 'relative': values in change are relative factors (e.g. 0.5 for 50% decrease)
%                                         2 or 'absolute': values are absolute values that overwrite the values in strct
% OUTPUT
% strct: updated neuron structure
% cell arrays of neuron structure are not supported so far
% Copyright by Marcel Beining <marcel.beining@gmail.com>

if nargin < 3 || isempty(mode) || (ischar(mode) && strcmpi(mode,'relative'))
    mode = 1;
else
    mode = 2;
end
changer = ones(1,2);
flag = 0;
if isfield(strct,'mech')
	ostrct = strct;
	strct = strct.mech;
	flag = 1;
end
if ~isempty(change) && isstruct(change) && isstruct(strct)
    
    chfield = fieldnames(change);
    field = fieldnames(strct);
    for f = 1:numel(field)
        if isstruct(strct.(field{f}))
            ffield = fieldnames(strct.(field{f}));
            inters = intersect(ffield,chfield);
            if ~isempty(inters)
                for in = 1:numel(inters)
                    changer(1) =  strct.(field{f}).(inters{in});
                    strct.(field{f}).(inters{in}) = changer(mode) * change.(inters{in});
                end
            end
            for ff = 1:numel(ffield)
                if isstruct(strct.(field{f}).(ffield{ff}))
                    inters = intersect(fieldnames(strct.(field{f}).(ffield{ff})),chfield);
                    if ~isempty(inters)
                        for in = 1:numel(inters)
                            changer(1) = strct.(field{f}).(ffield{ff}).(inters{in});
                            strct.(field{f}).(ffield{ff}).(inters{in}) = changer(mode) * change.(inters{in});
                        end
                    end                  
                end
            end
        end
    end
end
if flag
	ostrct.mech = strct;
	strct = ostrct;
end