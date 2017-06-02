function [bAPdisthm,mveloc_dend,mveloc_farax,mveloc_nearax,fig] = t2n_plotbAP(targetfolder_data,targetfolder_results,nneuron,ostruct)
if ~iscell(nneuron)
    nneuron = {nneuron};
end
thisdist = 185; % dist in µm where bAP is measured
if ~isfield(ostruct,'dist')
    ostruct.dist = 'Eucl.';
end
if ~isfield(ostruct,'relamp')
    ostruct.relamp = 0;
end
load(t2n_expcat(targetfolder_data,'Exp_bAP',nneuron{1}.experiment),'params')
if nargout == 0 || nargout > 4
    fig(1) = figure;clf,hold all
    if  exist(fullfile(params.path,'raw data','krueppel_data_fig_1d.dat'),'file')
        if ~ostruct.relamp
            data = importdata(fullfile(params.path,'raw data','krueppel_data_fig_1d.dat'));
            
        else
            data = importdata(fullfile(params.path,'raw data','krueppel_data_fig_1e.csv'));
        end
    else
        data = NaN(1,2);
    end
    if  exist(fullfile(params.path,'raw data','krueppel_data_fig_1f.csv'),'file')
        data2 = importdata(fullfile(params.path,'raw data','krueppel_data_fig_1f.csv'));
    else
        data2 = NaN(1,2);
    end
    if ostruct.relamp
        ylabel('Rel. amplitude')
        ylim([0 1.1])
        xlabel(sprintf('Rel %s distance to soma [µm]',ostruct.dist))
    else
        ylabel('Amplitude [mV]')
        ylim([0 120])
        xlabel(sprintf('%s distance to soma [µm]',ostruct.dist))
        
    end
    FontResizer
    
    fig(2) = figure;clf,hold all
    xlabel(sprintf('%s distance to soma [µm]',ostruct.dist))
    ylabel('Delay [ms]')
    FontResizer
    % bAP{t} 2nd dim: nodes_ind, time of max amp, PL to root, Eucl to root, max Voltage, baseline voltage before AP];
    for n = 1:numel(nneuron)
        fig(n+2) = figure;clf,hold all
        axis off
    end
end

for n = 1:numel(nneuron)
    load(t2n_expcat(targetfolder_data,'Exp_bAP',nneuron{n}.experiment),'bAP','plotvals','params','nodes','neuron','tree','tim')
    
    spiked = cellfun(@(x) any(x(:,5)>0),bAP);

    xlims = [Inf Inf];
    ylims = xlims;
    for t = 1:numel(tree)
        som = find(bAP{t}(:,1)==1);
        if isfield(tree{t},'col')
            col = tree{t}.col{1};
        else
            col = [1 0 0];
        end
        if strcmp(ostruct.dist,'PL.')
            L = Pvec_tree(tree{t});
        else
            L = eucl_tree(tree{t});
        end
        [~,spikeinitnode]=min(bAP{t}(:,2));
        %delay
        y = NaN(numel(tree{t}.Y),1);
        rax = find(strncmp(tree{t}.rnames,'axon',4));
        dendind = find(all(repmat(tree{t}.R,1,numel(rax))~=repmat(rax,numel(tree{t}.R),1),2));
        axind = find(any(repmat(tree{t}.R,1,numel(rax))==repmat(rax,numel(tree{t}.R),1),2));
        [dendind2,iad] = intersect(nodes{t},dendind);
        y(dendind2) = (bAP{t}(iad,2)-bAP{t}(som,2)); %time minus time at soma ( as in krueppel) %bAP{t}(:,1) == 1
        if nargout == 0 || nargout > 4
            figure(fig(2))
            if ostruct.relamp
                plotadjval(L/max(L(dendind2)),y,tree{t},col);
            else
                plotadjval(L,y,tree{t},col);
            end
            %bAP
            figure(fig(1))
            y = NaN(numel(tree{t}.Y),1);
            y(dendind2) = bAP{t}(iad,5)-bAP{t}(iad,6); %amplitude minus baseline
            if ostruct.relamp
                y = y/max(y);
                ipar = ipar_tree(tree{t});
                TP = T_tree(tree{t});
                ipar = ipar(TP,:);
                PL2 = NaN(numel(tree{t}.Y),1);
                for d = 2:numel(dendind2)
                    [x,~] = find(dendind2(d)==ipar);
                    PL2(dendind2(d)) = L(dendind2(d))/max(L(ipar(x,1)));
                end
                plotadjval(PL2,y,tree{t},col);
            else
                plotadjval(L,y,tree{t},col);
            end
            
            figure(fig(n+2));
            ptree = tran_tree(rot_tree(tran_tree(tree{t}),[],'-m3dY'),[350*t 300 0]);
            ptree.D(ptree.D<2) = 2;
            
            axind2 = intersect(nodes{t},axind);
            plotvals{t}(axind2) = NaN;
            xlims = [min(xlims(1),min(ptree.X(dendind))),max(xlims(2),max(ptree.X(dendind)))];
            ylims = [min(ylims(1),min(ptree.Y(dendind))),max(ylims(2),max(ptree.Y(dendind)))];
            plot_tree(ptree,plotvals{t});
        end
        
        % ind nodes, time of max amp, PL at nodes, eucl at nodes, max amplitude, baseline, time of halfmax amp
        ind = abs((bAP{t}(som,5)-bAP{t}(som,6))/2 - (bAP{t}(iad,5)-bAP{t}(iad,6))) <= 1; %index of half maximum
        
        if strcmp(ostruct.dist,'PL')
            bAPdisthm(t) = mean(bAP{t}(iad(ind),3)); % distance of half maximum amplitude
            ind = find(bAP{t}(:,3) - thisdist < 1);
            veloc = bAP{t}(:,3)./(bAP{t}(:,2)-bAP{t}(bAP{t}(:,1) == 1,2)); %L / Zeit die amp gebraucht hat von soma zu punkt (wie krueppel)
            veloc_dend = veloc(iad);
            veloc = abs(bAP{t}(:,3)-bAP{t}(spikeinitnode,3))./(bAP{t}(:,2+5)-bAP{t}(spikeinitnode,2+5)); %L / Zeit die amp gebraucht hat von spikeiniation zu punkt (wie kress)
        else
            bAPdisthm(t) = mean(bAP{t}(iad(ind),4)); % distance of half maximum amplitude
            ind = find(bAP{t}(:,4) - thisdist < 1);
            veloc = bAP{t}(:,4)./(bAP{t}(:,2)-bAP{t}(som,2)); %L / Zeit die amp gebraucht hat von soma zu punkt, hier maxamp als zeitpunkt  (wie krueppel)
            veloc_dend = veloc(iad);
            veloc = abs(bAP{t}(:,4)-bAP{t}(spikeinitnode,4))./(bAP{t}(:,2+5)-bAP{t}(spikeinitnode,2+5)); %L / Zeit die amp gebraucht hat von spikeiniation zu punkt, hier halfmax als Zeitpunkt (wie kress bzw SH08)
%             latency = (bAP{t}(:,2+5)-bAP{t}(som,2+5));
        end
        idpar = idpar_tree(tree{t});
        ind = intersect(iad,setdiff(ind,idpar(ind)));  % get only dendritic nodes and delete all direct parent nodes (due to rough distance search)
        bAPrelthisdist(t) = mean(bAP{t}(ind,5)-bAP{t}(ind,6))/(bAP{t}(som,5)-bAP{t}(som,6));
        bAPdelaythisdist(t) = mean(bAP{t}(ind,2))-bAP{t}(som,2);
        
        faraxind = axind(L(axind) > 100);
        [~,iaa] = intersect(nodes{t},faraxind);
        nearaxind = axind(L(axind) <= 100);
        [~,iaa2] = intersect(nodes{t},nearaxind);
        veloc_farax = veloc(iaa);
        veloc_nearax = veloc(iaa2);
        veloc_dend = veloc_dend(~isnan(veloc_dend) & ~isinf(veloc_dend)); % remove not measured and infinite values
        veloc_farax = veloc_farax(~isnan(veloc_farax) & ~isinf(veloc_farax)); % remove not measured and infinite values
        veloc_nearax = veloc_nearax(~isnan(veloc_nearax) & ~isinf(veloc_nearax)); % remove not measured and infinite values
        mveloc_dend(t) = mean(veloc_dend); % mean of velocity
        mveloc_farax(t) = mean(veloc_farax); % mean of velocity
        mveloc_nearax(t) = mean(veloc_nearax); % mean of velocity
    end
    if nargout == 0 || nargout > 4
        figure(fig(1))
        if ostruct.usemorph >= 4  % rat
            plot(data(:,1),data(:,2),'Marker','.','color','k','markersize',10,'linestyle','none') %*MRratioPL
        end
        xlim([0 400])
        FigureResizer(5,8)
        if ostruct.relamp
            ylim([0 1])
            tprint(fullfile(targetfolder_results,sprintf('bAP-rel-ampl_%s',nneuron{n}.experiment)),'-pdf')
        else
            ylim([0 150])
            tprint(fullfile(targetfolder_results,sprintf('bAP-ampl_%s',nneuron{n}.experiment)),'-pdf')
        end
        figure(fig(2))
        if ostruct.usemorph >= 4  % rat
            if ostruct.relamp
                plot(data2(:,1)/300,data2(:,2),'Marker','.','color','k','markersize',10,'linestyle','none') %*MRratioPL
            else
                plot(data2(:,1),data2(:,2),'Marker','.','color','k','markersize',10,'linestyle','none') %*MRratioPL
            end
        end
        ylim([-0.5 4.5])
        xlim([0 400])
        FigureResizer(5,8)
        
        tprint(fullfile(targetfolder_results,sprintf('bAP-del_%s',nneuron{n}.experiment)),'-pdf')
        figure(fig(n+2))
        ylim(ylims)
        xlim(xlims)
        ostruct.image = 1 ;
        FigureResizer(5,17,[],ostruct)
        c = colorbar;
        c.Limits =[-80,80];
        set(c,'Position',[0.93 0.35 0.02 0.4],'fontweight','bold','fontname','Arial')
        set(c,'YTick',[-80,0,80])
        tprint(fullfile(targetfolder_results,t2n_expcat('bAP-trees',nneuron{n}.experiment)),'-SHR-tif')
    end
    fprintf('Dendritic Velocity cell %d: %f µm/ms (time to max amp)\n',reshape(cat(1,(1:numel(mveloc_dend)),mveloc_dend),1,numel(mveloc_dend)*2))
    fprintf('Far Axonal Velocity cell %d: %f µm/ms (time to half-max amp)\n',reshape(cat(1,(1:numel(mveloc_farax)),mveloc_farax),1,numel(mveloc_farax)*2))
    fprintf('Near Axonal Velocity cell %d: %f µm/ms (time to half-max amp)\n',reshape(cat(1,(1:numel(mveloc_nearax)),mveloc_nearax),1,numel(mveloc_nearax)*2))
    % hide p, set colorbar to border and save
    
    if ~all(spiked)
        display('CAUTION: Not all cells spiked!')
    end
    fprintf('Mean voltage attenuation @ %d µm: %g +- %g %% (s.e.m.)\n',thisdist,mean(bAPrelthisdist)*100,std(bAPrelthisdist)/sqrt(numel(tree))*100)
    if (nargout == 0 || nargout > 4) && ~isnan(data(1))
        fprintf('Mean voltage attenuation in exp @ 185 µm: %g +- %g %% (s.e.m.)\n',mean(data(17:20,2))/mean(cellfun(@(x) x(1,5)-x(1,6),bAP))*100,std(data(17:20,2))/mean(cellfun(@(x) x(1,5)-x(1,6),bAP))/sqrt(4)*100)  % bAP at 185 µm in exp
    end
    fprintf('Mean delay @ %d µm: %g +- %g ms (s.e.m.)\n',thisdist,mean(bAPdelaythisdist),std(bAPdelaythisdist)/sqrt(numel(tree)))
    
end

function h = plotadjval(x,y,tree,col)

tree.Z(:) = 0;
tree.X = x;
tree.Y = y;
if sum(isnan(tree.Y)) ~= 0
    tree = delete_tree(tree,find(isnan(tree.Y)));
end
treecol = repmat(col,numel(tree.Y),1);
h = plot_tree(tree,treecol,[],[],[],'-2l');
axis normal