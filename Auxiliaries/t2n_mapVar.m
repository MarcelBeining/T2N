function [varVec] = t2n_mapVar(tree,neuron,var,show)
% This function maps the values of a given NEURON variable onto each
% morpholgy and optionally returns these values, too;
% INPUT
% tree      TREES toolbox morphologies
% neuron    t2n neuron structure or cell array of structure (see
%           documentation)
% var       string with valid NEURON variable name, such as 'g_pas'
% show      (optional) boolean if morphologies should be plotted (default)
%
% OUTPUT
% varVec    m-by-n cell array with mapped values for m neuron
%           specifications and n trees

if nargin < 4 || isempty(show)
    show = 1;
end
if isstruct(tree)
    tree = {tree};
end
if isstruct(neuron)
    neuron = {neuron};
end
varSplit = regexp(var,'_','split');  % split variable into mechanism name and variable name
varVec = cell(numel(neuron),numel(tree));
for n = 1:numel(neuron) % go through all neuron definitions
    for t = 1:numel(tree)  % go through all morphologies
        varVec{n,t} = NaN(numel(tree{t}.X),1);  % initialize the vector as nans
        if isfield(neuron{n},'mech')   % check for mechanism definition
            regions = fieldnames(neuron{n}.mech{t}); % get region names
            if any(strcmp(regions,'all'))   % variable was set in all nodes
                if isfield(neuron{n}.mech{t}.all,varSplit{2})  % check if mechanism exists at that region
                    if isfield(neuron{n}.mech{t}.all.(varSplit{2}),varSplit{1}) % check if variable is defined at that location
                        varVec{n,t}(:) = neuron{n}.mech{t}.all.(varSplit{2}).(varSplit{1}); % save value in all nodes
                    else
                        varVec{n,t} = Inf;  % mechanism seems to be defined but not the variable. remember that
                    end
                end
               regions = regions(~strcmp(regions,'all')); % delete region 'all' from regions
            end
            regions = intersect(regions,tree{t}.rnames);  % remove all nondefined regions
            for r = 1:numel(regions)  % go through regions
                if isfield(neuron{n}.mech{t}.(regions{r}),varSplit{2})  % check if mechanism exists at that region
                    if isfield(neuron{n}.mech{t}.(regions{r}).(varSplit{2}),varSplit{1}) % check if variable is defined at that location
                        ind = tree{t}.R == find(strcmp(tree{t}.rnames,regions{r}));  % get index to all nodes that belong to this region
                        varVec{n,t}(ind) = neuron{n}.mech{t}.(regions{r}).(varSplit{2}).(varSplit{1}); % save value in these nodes
                    else
                        varVec{n,t} = Inf;  % mechanism seems to be defined but not the variable. remember that
                    end
                end
            end
            if isfield(neuron{n}.mech{t},'range')  % check for a range structure
                mechanisms = intersect(varSplit{2},fieldnames(neuron{n}.mech{t}.range)); % check if our mechanism exists in the range variable
                if ~isempty(mechanisms)
                    vars = intersect(varSplit{1},fieldnames(neuron{n}.mech{t}.range.(mechanisms))); % check if our variable exists in the range variable
                    if ~isempty(vars)
                        ind = ~isnan(neuron{n}.mech{t}.range.(mechanisms).(vars));  % get only the not nan values (nan means use the default value or value defined by the region spec)
                        varVec{n,t}(ind) = neuron{n}.mech{t}.range.(mechanisms).(vars)(ind); % save the not nan values
                    end
                end
            end
        end
        if any(isinf(varVec{n,t}))
            warning('Caution! Values of %s seem to have not been defined at some nodes! NEURON uses the mechanisms default value (which T2N does not know), hence these values cannot be displayed.',var)
            varVec{n,t}(isinf(varVec{n,t})) = NaN;
        end
        if show  % map values on tree
            figure
            if isfield(tree{t},'name')
                tname = tree{t}.name;
            else
                tname = sprintf('tree %d',t);
            end
            title(sprintf('Neuron Simulation: %d, tree: "%s", variable "%s"',n,tname,strrep(var,'_','\_')))
            plot_tree(tree{t},varVec{n,t})
            lims = [min(varVec{n,t}),max(varVec{n,t})];
            if diff(lims) == 0  % limits HAVE to be different from each other
                lims = [lims(1)-1e-9,lims(1)+1e-9];
            end
            set(gca,'CLim',lims)
        end
    end
end

end

