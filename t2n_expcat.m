function Filename = t2n_expcat(varargin)
% concatenates strings to experiment name
if strcmp(varargin{end-1}(1),'_')
    Filename = strcat(varargin{end-1},varargin{end}) ;
else
    Filename = strcat(varargin{end-1},'_',varargin{end});
end

Path = '';
for n = 1:numel(varargin)-2
    Path = strcat(Path,'/',varargin{n});
end

Filename = strcat(Path(2:end),'/',Filename,'.mat');
