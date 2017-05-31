function [ vec, tvecout ] = t2n_get(out,par,arg,typ)
% get a parameter at a specific time point/period
% INPUT :   out     is the output structure of t2n
%           par     is the parameter to obtain (string)
%           arg     is either a time point [ms] or a start and end time
%                   point [ms] or the string name of a function that should be
%                   applied on the recorded values at each node
%           typ     specifies if the parameter is from the cell ('cell',
%                   default) or a point process (e.g. 'IClamp')
% currently not supporting local_dt

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
        vec{o}{t} = NaN(numel(out{o}.record{t}.(typ).(par)),numl);
        for n = 1:numel(out{o}.record{t}.(typ).(par))
            vec{o}{t}(n,:) = fHandle(out{o}.record{t}.(typ).(par){n}) ;
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

