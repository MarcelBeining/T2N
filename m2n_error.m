function out = m2n_error(out,outoptions,errorcode)
if nargin < 3 || isempty(errorcode)
    errorcode = 1;
end
if nargin < 2 || isempty(outoptions)
    outoptions.nocell = 0;
end

for n = 1:numel(out)
        out{n}.error = errorcode;
end

if outoptions.nocell  % make output a structure not cell since only one simulation was calculated..
   out = out{1}; 
end