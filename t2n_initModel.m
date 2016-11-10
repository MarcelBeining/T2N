function [tree,params,neuron,treename] = t2n_initModel(params,ostruct)

if ~isfield( ostruct,'reducecells')
    ostruct.reducecells = 0;
end
if ~isfield( ostruct,'channelblock')
    ostruct.channelblock = {};
end
if ~isfield( ostruct,'newborn')
    ostruct.newborn = 0;
end
if ~isfield( ostruct,'forcecalcload')
    ostruct.forcecalcload = 0;
end



origfolder = pwd;

if ostruct.vmodel == 0
    mechoptions = '-p';
elseif ostruct.vmodel > 0
    mechoptions = '-a-p-zz';
else
    mechoptions = '-o-a-p';
end
if ostruct.newborn
    mechoptions = strcat(mechoptions,'-y');
end
if ostruct.ratadjust
    mechoptions = strcat(mechoptions,'-ra');
end

cd(params.path)
if ostruct.ratadjust
    str = '_ratadjust';
else
    str = '';
end
if ~isnan(ostruct.usemorph)
    
    switch ostruct.usemorph
        case 1    % SH07 cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed.mtr'));
            end
            
            params.tname = 'SH07all2';
        case 2   % synth adult mouse cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon.mtr'));
            end
            params.tname = 'mouse_matGC_art';
            params.exchfolder = 't2nexchange_aGCmorphsim';
            
        case 3  % synth young mouse cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon.mtr'));
            end
            params.tname = 'mouse_abGC_art';
            params.exchfolder = 't2nexchange_aGCmorphsim2';
        case 4  %  adult rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,sprintf('/morphos/Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/Beining_AAV_contra_MLyzed_axon.mtr'));
            end
            params.tname = 'rat_mGC_Beining';
            params.exchfolder = 't2nexchange_aGCmorphsim6';
        case 5 % synth adult rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,sprintf('/morphos/rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/rat_AAVart_old_pruned_axon.mtr'));
            end
            params.tname = 'rat_matGC_art';
            params.exchfolder = 't2nexchange_aGCmorphsim3';
        case 6  % synth young rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,sprintf('/morphos/rat_RVart_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/rat_RVart_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/rat_RVart_pruned_axon.mtr'));
            end
            params.tname = 'rat_abGC_art';
            params.exchfolder = 't2nexchange_aGCmorphsim4';
        case 7  % Claiborne rat cells
            if ostruct.adjustloads && ~ostruct.forcecalcload && exist(fullfile2(params.path,sprintf('/morphos/Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/Claiborne_male_MLyzed_axon.mtr'));
            end
            params.tname = 'rat_mGC_Claiborne';
            params.exchfolder = 't2nexchange_aGCmorphsim5';

    end
    neuron.experiment = params.tname;
    
        ntree = min(numel(tree),15);
        tree=tree(1:ntree);
        if ostruct.usecol
            if any(ostruct.usemorph == [1,4]) % reconstructed morphologies
                colors = colorme(ntree,'-grb');
            else
                colors = colorme(ntree,'-grg');
            end
        else
            colors = colorme(ntree);
        end
        for t = 1:numel(tree)
            tree{t}.col = colors(t);
        end
else
    [tree,treename,treepath]=load_tree;
    if isempty(tree)
        return
    end

    colors = colorme(numel(tree));
    for t = 1:numel(tree)
        tree{t}.col = colors(t);
    end
    params.tname = treename(1:end-4);
end


if ~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(params.morphfolder,[x.NID,'.hoc']),'file'),tree))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        tree = t2n_writetrees(params,tree,[],fullfile(treepath,treename));
    end
end
if isempty(tree)
    return
end
if isstruct(tree)
    tree={tree};
elseif iscell(tree{1})
    tree=tree{1};
end

cd(params.path)

for t = 1:numel(tree)
    tree{t} = sort_tree(tree{t},'-LO');
    
    neuron.mech{t} = [];
    if ostruct.scalespines
        neuron.mech{t} =  cat_struct(AH_active(tree{t},mechoptions),spinedens_mGC(ostruct.scalespines*0.9));
    else
        neuron.mech{t} =  cat_struct(AH_active(tree{t},mechoptions));
    end
    if ~isempty(strfind(mechoptions,'-Ca')) && isempty(strfind(mechoptions,'-Kv7')) ||  ~isempty(strfind(mechoptions,'buff'))
        if ~isempty(strfind(mechoptions,'-Nabuff'))
            neuron = setionconcentrations(neuron,'Mongiat3');
        else
            neuron = setionconcentrations(neuron,'Mongiat');
        end
    else
        neuron = setionconcentrations(neuron,'Mongiat');
    end
    if ~isfield(tree{t},'col')
        tree{t}.col{1} = rand(1,3);
    end
    if ostruct.noise ~= 0
        neuron.pp{t}.InGauss = struct('node',1,'mean',0.01,'stdev',0.01,'del',0,'dur',1e9);
    end
end

if ostruct.scalespines
    neuron = scale_spines(neuron);
end
%%% This is the Hay 2013 implementation of adjusting soma and AIS
%%% conductance according to dendritic morphology
if ostruct.adjustloads
    if any(cellfun(@(x) ~isfield(x,'Rho_soma') | ~isfield(x,'Rho_AIS'),tree)) || ostruct.forcecalcload
        tree = calculate_loads(params,neuron,tree);
        save_tree(tree,fullfile(treepath,[treename(1:end-4),sprintf('_loadadjusted%s.mtr',str)]));
    end
    if ostruct.usemorph >= 4      %rat
        neuron = adjust_loads(neuron,tree,'r',ostruct);
    else        %mouse
        neuron = adjust_loads(neuron,tree,'m');
    end
end

if ostruct.reducecells
    if ostruct.usemorph ~= 1
        tree=tree((1:3)+2);
        neuron.mech = neuron.mech((1:3)+2);
    else
        tree=tree(3);
        neuron.mech = neuron.mech(3);
    end
    neuron.experiment = strcat(neuron.experiment,'_reduceNs');
end
if ostruct.ratadjust
    neuron.experiment = strcat(neuron.experiment,'_ratadjust');
end

if ~isempty(ostruct.channelblock)
    neuron = blockchannel(neuron,ostruct.channelblock,ostruct.blockamount,ostruct.specify);
    display(sprintf('Channel(s) %s blocked!',cell2mat(ostruct.channelblock)))
end

treename = fullfile(treepath,treename);
cd(origfolder)
display(sprintf('Model initialized...%s',neuron.experiment))

