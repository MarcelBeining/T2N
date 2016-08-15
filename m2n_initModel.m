function [tree,params,neuron] = m2n_initModel(params,ostruct)

if ~isfield( ostruct,'reducecells')
    ostruct.reducecells = 0;
end
if ~isfield( ostruct,'channelblock')
    ostruct.channelblock = {};
end
if ~isfield( ostruct,'newborn')
    ostruct.newborn = 0;
end
origfolder = pwd;

if ostruct.vmodel == 0
    mechoptions = '-p-g';
elseif ostruct.vmodel == 1
    mechoptions = '-a-p-n-g-8st-sh3-Ca';%-sh
elseif ostruct.vmodel == 2
    mechoptions = '-a-p-n2-g-8st-sh3-Ca';%-sh
elseif ostruct.vmodel == 3
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2';%-sh
elseif ostruct.vmodel == 4
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2-TWIK';%-sh
elseif ostruct.vmodel == 5
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2-TWIK-NaKbuff';%-sh
elseif ostruct.vmodel == 6
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2-TWIK-Kv7';%-sh
elseif ostruct.vmodel == 7
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2-TWIK-Kv7-nl';%-sh
elseif ostruct.vmodel == 8
    mechoptions = '-a-p-n2-g-8st-sh3-Ca2-TWIK-Kv7-nl';%-sh %same  but not in spinedens
    
else
    mechoptions = '-o-a-p';
    Mname = 'AH99';
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
if ~isnan(ostruct.usemorphmodel)
    
    switch ostruct.usemorphmodel
        case 1
            if ostruct.adjustloads && exist(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_AAVart_old_pruned_axon.mtr'));
            end
            params.tname = 'mouse_matGC_art';
            params.exchfolder = 'm2nexchange_aGCmorphsim';
            
        case 2
            if ostruct.adjustloads && exist(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/mouse_RVart_pruned_axon.mtr'));
            end
            params.tname = 'mouse_abGC_art';
            params.exchfolder = 'm2nexchange_aGCmorphsim2';
        case 3
            if ostruct.adjustloads && exist(fullfile2(params.path,sprintf('/morphos/rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/rat_AAVart_old_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/rat_AAVart_old_pruned_axon.mtr'));
            end
            params.tname = 'rat_matGC_art';
            params.exchfolder = 'm2nexchange_aGCmorphsim3';
        case 4
            if ostruct.adjustloads && exist(fullfile2(params.path,sprintf('/morphos/rat_RVart_pruned_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/rat_RVart_pruned_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/rat_RVart_pruned_axon.mtr'));
            end
            params.tname = 'rat_abGC_art';
            params.exchfolder = 'm2nexchange_aGCmorphsim4';
        case 5
            if ostruct.adjustloads && exist(fullfile2(params.path,sprintf('/morphos/Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/Claiborne_male_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/Claiborne_male_MLyzed_axon.mtr'));
            end
            params.tname = 'rat_mGC_Claiborne';
            params.exchfolder = 'm2nexchange_aGCmorphsim5';
        case 6
            if ostruct.adjustloads && exist(fullfile2(params.path,sprintf('/morphos/Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,sprintf('/morphos/Beining_AAV_contra_MLyzed_axon_loadadjusted%s.mtr',str)));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/Beining_AAV_contra_MLyzed_axon.mtr'));
            end
            params.tname = 'rat_mGC_Beining';
            params.exchfolder = 'm2nexchange_aGCmorphsim6';
        case 0
            if ostruct.adjustloads && exist(fullfile2(params.path,'/morphos/SH_07_MLyzed_new3_soma_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_MLyzed_new3_soma_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_MLyzed_new3_soma.mtr'));
            end
            
            params.tname = 'SH07';
        case 0.5
            if ostruct.adjustloads && exist(fullfile2(params.path,'/morphos/SH_07_all_repairedandsoma_MLyzed_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsoma_MLyzed_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsoma_MLyzed.mtr'));
            end
            
            params.tname = 'SH07all';
            %             params.exchfolder = 'm2nexchange_aGC_SH07';
        case 0.6
            if ostruct.adjustloads && exist(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'),'file')
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed_loadadjusted.mtr'));
            else
                [tree,treename,treepath]=load_tree(fullfile2(params.path,'/morphos/SH_07_all_repairedandsomaAIS_MLyzed.mtr'));
            end
            
            params.tname = 'SH07all2';
            %             params.exchfolder = 'm2nexchange_aGC_SH07';
    end
    neuron.experiment = params.tname;
    if ostruct.usemorphmodel > 0
        ntree = min(numel(tree),15);
        tree=tree(1:ntree);
        %         Mname = strcat(Mname,sprintf('_morphmodel%d',options.usemorphmodel));
        if ostruct.usecol
            if ostruct.usemorphmodel < 1 || ostruct.usemorphmodel > 5 % reconstructed morphologies
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
        colors = {'dim blue','blue','light blue'};
        for t = 1:numel(tree)
            tree{t}.col = colorme(colors{t});
        end
    end
else
    %     if ~exist('treename','var') || (exist('treename','var') && isempty(strfind(treename,'SH_07_MLyzed'))) || isempty(tree) || numel(tree) == 1
    [tree,treename,treepath]=load_tree;
    if isempty(tree)
        return
    end
    %         if ~isempty(strfind(treename,'new2'))
    % %             tree = tree(cellfun(@(x) any(strcmp(usetrees,x.name)),tree));
    %
    %         end
    colors = colorme(numel(tree));
    for t = 1:numel(tree)
        tree{t}.col = colors(t);
    end
    params.tname = treename(1:end-4);
end


if ~all(cellfun(@(x) isfield(x,'NID'),tree)) || ~all(cellfun(@(x) exist(fullfile(params.morphfolder,[x.NID,'.hoc']),'file'),tree))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        tree = m2n_writetrees(params,tree,[],fullfile(treepath,treename));
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
% neuron.tree = 1:numel(tree);
for t = 1:numel(tree)
    tree{t} = sort_tree(tree{t},'-LO');
    
    %     tree{t}.col = colorgb(t);
    neuron.mech{t} = [];
    %         [neuron.mech{t}, addsurf(t)] =  spinedens_mGC([],tree{t});
    %         tsurf(t) = sum(surf_tree(tree{t})) / 10^8;     % cm²
    if ostruct.scalespines
        neuron.mech{t} =  cat_struct(AH_active(tree{t},mechoptions,t),spinedens_mGC(ostruct.vmodel));%,SH_na8st(tree{t}));  %mature spine density and AH model of active cell !!!!!!!!!!!!!!!!!!!
    else
        neuron.mech{t} =  cat_struct(AH_active(tree{t},mechoptions,t));%,SH_na8st(tree{t}));  %mature spine density and AH model of active cell !!!!!!!!!!!!!!!!!!!
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
    %     if options.newborn
    %        neuron.mech{t}.all.k_ion.ek = -80;
    %     end
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
    if ostruct.usemorphmodel >= 3      %rat
        neuron = adjust_loads(neuron,tree,'r',ostruct);
    else        %mouse
        neuron = adjust_loads(neuron,tree,'m');
    end
end

if ostruct.reducecells
    if ostruct.usemorphmodel >= 0.1% || isempty(strfind(params.tname,'SH07'))
        tree=tree(1:3);
        neuron.mech = neuron.mech(1:3);
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
    neuron = blockchannel(neuron,ostruct.channelblock,ostruct.blockamount);
    display(sprintf('Channel(s) %s blocked!',cell2mat(ostruct.channelblock)))
end
if ostruct.cAMP ~= 0
    neuron = changeoptions.cAMP(neuron,ostruct.cAMP);
end

cd(origfolder)
display(sprintf('Model initialized...%s',neuron.experiment))

