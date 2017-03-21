function strct = cat_struct(varargin)
% 
% this function concatenates two structures which have the same fields iteratively
% this function is part of the T2N package
% INPUT
% several structures with same field names
% OUTPUT
% concatenated structure
%
% Copyright by Marcel Beining <marcel.beining@gmail.com>

strct = struct();
for v = 1:nargin
    fields = fieldnames(varargin{v});
    for f = 1:numel(fields)
        if isfield(strct,fields{f})
            strct.(fields{f}) = cat_struct(strct.(fields{f}),varargin{v}.(fields{f}));
        else
            strct.(fields{f}) = varargin{v}.(fields{f});
        end
    end
end