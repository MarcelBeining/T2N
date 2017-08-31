function [ vec, tvecout ] = t2n_get(out,par,arg,typ)
% This function can be used to obtain a recorded parameter during a specific 
% time point/period of the simulation.
%
% INPUTS
% out     the output structure of the t2n main function
% par     string of the the parameter to obtain
% arg     can be either
%         - a scalar defining time point [ms] at which the parameter should be returned
%         - a 1x2 vector with the start and end time point [ms] between which the parameter should be returned
%         - the string name of a function that should be applied on the recorded parameter values at each node
% typ     string that specifies if the recorded parameter is from the cell ('cell',default) or from a point process (e.g. 'IClamp')
%
% OUTPUTS
% vec       cell array or vector returning the desired parameter at simulation (first level) and node (second level)  
% tvecout   time vector [ms] if arg was no string
%
% NOTE
% This function is currently not supporting NEURON's local_dt feature
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 4 || isempty(typ)
    typ = 'cell';
end
nocellflag = 0;
if ~iscell(out)
    out = {out};
    nocellflag = 1;
end
vec = cell(numel(out),1);
tvecout = vec;
for o = 1:numel(out)
    switch class(arg)
        case 'char'
            fHandle = str2func(arg);
            numl = 1;
        otherwise
            switch numel(arg)
                case 1
                    tvec = find(out{o}.t >= arg,1,'first');
                    numl = 1;
                case 2
                    tvec = out{o}.t >= arg(1) & out{o}.t <= arg(2);
                    numl = sum(tvec);
            end
            fHandle = @(x) x(tvec);
    end
    
    %     vec{p} = NaN(numel(out.record),max(cellfun(@(x) numel(x.(typ).(par{p})),out.record)))  % initialize the matrix with N
    for t = 1:numel(out{o}.record)
        if ~isempty(out{o}.record{t})
            vec{o}{t} = NaN(numel(out{o}.record{t}.(typ).(par)),numl);
            for n = 1:numel(out{o}.record{t}.(typ).(par))
                vec{o}{t}(n,:) = fHandle(out{o}.record{t}.(typ).(par){n}) ;
            end
        end
    end
    if ~isa(arg,'char')
        tvecout{o} = out{o}.t(tvec);
    else
        tvecout{o} = NaN;
    end
end

if nocellflag || numel(out) == 1
    vec = vec{1};
    tvecout = tvecout{1};
    if all(cellfun(@numel,vec) == 1)
        vec = cell2mat(vec);
    end
elseif all(cellfun(@(x) all(cellfun(@numel,x) == 1),vec))
    vec = cellfun(@cell2mat,vec,'UniformOutput',0);
end


end

