function strct = t2n_catStruct(varargin)
% This function concatenates two structures which have the same fields.
%
% INPUTS
% varargin          several structures with same field names
%
% OUTPUTs
% strct             concatenated structure
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

strct = struct();
for v = 1:nargin
    fields = fieldnames(varargin{v});
    for f = 1:numel(fields)
        if isfield(strct,fields{f})
            strct.(fields{f}) = t2n_catStruct(strct.(fields{f}),varargin{v}.(fields{f}));
        else
            strct.(fields{f}) = varargin{v}.(fields{f});
        end
    end
end