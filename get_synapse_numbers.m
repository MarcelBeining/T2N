function [nsyn,synids] = get_synapse_numbers(tree,syn_dens,regions)

if nargin < 3
    regions = tree.rnames;
end
if ~iscell(regions)
    regions = {regions};
end
if nargin < 2
    syn_dens = 1;
end
if numel(syn_dens) == 1 && numel(regions) > 1
    syn_dens = repmat(syn_dens,numel(regions),1);
end

nsyn = zeros(numel(tree.X),1);
synids = [];
for r = 1:numel(regions)
    
    ind = tree.R==find(strcmp(tree.rnames,regions{r}));
    
    sect = dissect_tree(tree);
    sect = sect(ind(sect(:,2)),:);
    ipar = ipar_tree(tree,[],sect(:,2));
    for n = 1:size(ipar,1)
        ipar(n,find(ipar(n,:)==sect(n,1)):end) = 1;  % make branch point and following nodes to the root refer to node 1 (root)
    end
    ipar = ipar(:,any(ipar ~= 1,1));  % delete unneeded columns
    len = len_tree(tree);        % get lengths of node segments
    synm = len(ipar)*syn_dens;   % calculate matrix with number of synapses per node
    synm(synm==0) = NaN;
    synm = round(synm*10)/10;    % round number of synapses
    flag = false;   % check flag if cumsum has been performed
    while 1 %any(synm(:)>=1)
        synids = cat(1,synids,ipar(synm>=1));
        nsyn(ipar(synm>=1)) = nsyn(ipar(synm>=1)) + 1;
        synm(synm>=1) = synm(synm>=1) -1;
        if ~any(synm(:)>=1)   % this has to be done to make it able to have non-integer spine densities
            if ~flag  
                for n = 2:size(synm,2)  % go through syn matrix and add up synapse numbers each time until value exceeds 1
                    add = synm(:,n-1);
                    add(add >= 1) = add(add >= 1) -1;
                    synm(:,n) = round((synm(:,n) + add)*10)/10;
                end
                flag = true;
            else
                break
            end
        end
    end
end