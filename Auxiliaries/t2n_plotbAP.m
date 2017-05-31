function [bAPdisthm,mveloc_dend,mveloc_farax,mveloc_nearax,fig] = t2n_plotbAP(targetfolder_data,targetfolder_results,nneuron,ostruct)
style = {'-','--',':','-.'};
if ~iscell(nneuron)
    nneuron = {nneuron};
end
thisdist = 185; % dist in µm where bAP is measured
if ~isfield(ostruct,'show')
    ostruct.show = 0;
end
if ~isfield(ostruct,'dist')
    ostruct.dist = 'Eucl.';
end
if ~isfield(ostruct,'relamp')
    ostruct.relamp = 0;
end
load(t2n_expcat(targetfolder_data,'Exp_bAP',nneuron{1}.experiment),'params')
if ostruct.show 
    if any(ostruct.show ==4)
        fig(8) = figure;hold on
        fig(6) = figure;hold on
    end
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
%     if  exist(fullfile(params.path,'raw data','AxonalVelocity_ratGC.csv'),'file')
%         data3 = importdata(fullfile(params.path,'raw data','AxonalVelocity_ratGC.csv'));
%     else
%         data3 = NaN(1,2);
%     end
    % title('bAP vs Plen')
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
%     fig(3) = figure;clf,hold all
%     xlabel(sprintf('%s distance to soma [µm]',ostruct.dist))
%     ylabel('Axonal Conduction Velocity [m/s]')
%     % line([0 1000],[149.3 149.3]/1000,'LineStyle','--','Color',[0.5 0.5 0.5])
%     plot(data3(:,1),data3(:,2),'Marker','.','color','k','markersize',10,'linestyle','none')
    % bAP{t} 2nd dim: nodes_ind, time of max amp, PL to root, Eucl to root, max Voltage, baseline voltage before AP];
    for n = 1:numel(nneuron)
        fig(n+3) = figure;clf,hold all
        axis off
    end
end

for n = 1:numel(nneuron)
    load(t2n_expcat(targetfolder_data,'Exp_bAP',nneuron{n}.experiment),'bAP','plotvals','params','nodes','neuron','tree','CaNodes','mCai','stdCai','mV','stdV','tim','tw','maxcai')
    
    spiked = cellfun(@(x) any(x>0),mV(:,1));
    
    
    
    
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
        if ostruct.show
            figure(fig(2))
            if ostruct.relamp
                plotadjval(L/max(L(dendind2)),y,tree{t},col);
            else
                plotadjval(L,y,tree{t},col);
            end
            %bAP
            figure(fig(1))
            %     plot(bAP{t}(:,3),bAP{t}(:,5)-bAP{t}(:,6),'Marker','.','Color',tree{t}.col{1},'markersize',4,'linestyle','none','linewidth',1)
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
            
            figure(fig(n+3));
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
        
%         if ostruct.show
%             figure(fig(3));
%             y = NaN(numel(tree{t}.Y),1);
%             y(axind2) = veloc(iaa)/1000;
%             plotadjval(L,y,tree{t},col);
%         end
        %     plot(PL(iaa),veloc(iaa),'x')
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
        
        
        if any(ostruct.show == 4)
            for f =1:numel(CaNodes{t})
                if 1%t == 1
                    figure(fig(6))
                    hold on
                    if ~all(stdCai{t,f}==0)
                        hp = patch ([tim', (fliplr (tim'))], [(mCai{t,f} + stdCai{t,f}), (fliplr (mCai{t,f} - stdCai{t,f}))],[0 0 0]);
                        
                        set (hp, 'facealpha', 0.2, 'edgecolor', 'none','facecolor',tree{t}.col{1});%col{f})
                    end
                    plot(tim,mCai{t,f},'Color',tree{t}.col{1},'LineStyle',style{f});%)
                    figure(fig(8))
                    hold on
                    if ~all(stdV{t,f}==0)
                        hp = patch ([tim' (fliplr (tim'))], [(mV{t,f} + stdV{t,f}) (fliplr (mV{t,f} - stdV{t,f}))],[0 0 0]);
                        
                        set (hp, 'facealpha', 0.2, 'edgecolor', 'none','facecolor',tree{t}.col{1})
                    end
                    plot(tim,mV{t,f},'Color',tree{t}.col{1},'LineStyle',style{f})
                    
                end
            end
        end
    end
    if ostruct.show
        figure(fig(1))
        if ostruct.usemorph >= 4  % rat
            plot(data(:,1),data(:,2),'Marker','.','color','k','markersize',10,'linestyle','none') %*MRratioPL
        end
        xlim([0 400])
        FigureResizer(5,8)
        if ostruct.relamp
            ylim([0 1])
            tprint(fullfile(targetfolder_results,sprintf('Fig.2-bAP-rel-ampl_%s',nneuron{n}.experiment)),'-pdf')
        else
            ylim([0 150])
            tprint(fullfile(targetfolder_results,sprintf('Fig.2-bAP-ampl_%s',nneuron{n}.experiment)),'-pdf')
        end
        
        %     tprint(fullfile(targetfolder_results,'Fig.2-bAP-ampl'),'-png')
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
        
        tprint(fullfile(targetfolder_results,sprintf('Fig.2-bAP-del_%s',nneuron{n}.experiment)),'-pdf')
        %     tprint(fullfile(targetfolder_results,'Fig.2-bAP-del'),'-png')
        figure(fig(n+3))
        ylim(ylims)
        xlim(xlims)
        %     ylim([250 500])
        %     xlim([200 1200])
        ostruct.image = 1 ;
        FigureResizer(5,17,[],ostruct)
        c = colorbar;
        c.Limits =[-80,80];
%         pos = get(c,'Position');
        set(c,'Position',[0.93 0.35 0.02 0.4],'fontweight','bold','fontname','Arial')
%         yl = get(c,'YLim');
        set(c,'YTick',[-80,0,80])
%         set(c,'YTick',[ceil(yl(1)),0,floor(yl(2))])
        tprint(fullfile(targetfolder_results,t2n_expcat('Fig.2-bAP-trees',nneuron{n}.experiment)),'-SHR-tif')
        %     tprint(fullfile(targetfolder_results,'Fig.2-bAP-trees'),'-SHR-png')
%         figure(fig(3))
%         ylim([0 0.5])
%         xlim([0 400])
%         FigureResizer(5,8)
%         tprint(fullfile(targetfolder_results,t2n_expcat('Fig.2-bAP-axon',nneuron{n}.experiment)),'-SHR-pdf')
    end
    %     set(p,'Visible','off')
    fprintf('Dendritic Velocity cell %d: %f µm/ms (time to max amp)\n',reshape(cat(1,(1:numel(mveloc_dend)),mveloc_dend),1,numel(mveloc_dend)*2))
    fprintf('Far Axonal Velocity cell %d: %f µm/ms (time to half-max amp)\n',reshape(cat(1,(1:numel(mveloc_farax)),mveloc_farax),1,numel(mveloc_farax)*2))
    fprintf('Near Axonal Velocity cell %d: %f µm/ms (time to half-max amp)\n',reshape(cat(1,(1:numel(mveloc_nearax)),mveloc_nearax),1,numel(mveloc_nearax)*2))
    %     close(778:780)
    % hide p, set colorbar to border and save
    
    if ostruct.show
        if any(ostruct.show==4)
            figure(fig(6))
            %     subplot(2,1,1)
            %     legend([p2{1}(end),p2{2}(end),p2{3}(end)],'Soma','Proximal','Distal')
            ylabel('Ca concentration [nM]')
            xlim([0 200])
            xlabel('Time [ms]')
            ylim([0, 1000]);%ceil(max(max(cat(2,out.record{t}.cell.cai{CaNodes{t}{2}})))*1e6)])
            figure(fig(8))
            %     subplot(2,1,2)
            xlim([30 60])
            ylim([-80 50])
            ylabel('Membrane Potential [mV]')
            xlabel('Time [ms]')
            
            %     tprint(fullfile(targetfolder_results,t2n_expcat(sprintf('Fig.X-CaDyn-tree%d',t),neuron{n}.experiment)),'-HR-png');
            
            figure(fig(6))
            FontResizer
            FigureResizer(5,8)
            %     tprint(fullfile(targetfolder_results,t2n_expcat('Fig.4-CaDyn',neuron{n}.experiment)),'-pdf');
            figure(fig(8))
            FontResizer
            FigureResizer(5,8)
            %     tprint(fullfile(targetfolder_results,t2n_expcat('Fig.4-CaDynV',neuron{n}.experiment)),'-pdf');
        end
        gstruct.ugroupdef = {{'Data Prox.','Proximal','Distal','Soma','Axon','Data MFB'}};
        if ~isempty(strfind(nneuron{n}.experiment,'_art'))
            gstruct.col = colorme('Black','Dim green','Dim green','Dim green','Dim green','Black');
        else
            gstruct.col = colorme('Black','Dim blue','Dim blue','Dim blue','Dim blue','Black');
        end
        gstruct.group{1} = 1:6;
        data_stocca = [0.23*1000,nanmean(nanmean(tw,3),1),43];
        sem = [0.03*1000,nanstd(nanmean(tw,3),[],1)/sqrt(numel(tree)),6];
        if ostruct.usemorph >= 4 && ~ostruct.newborn  % rat mature..add data
            data_stocca = data_stocca([1,4,5,3,2,6]);
            sem = sem([1,4,5,3,2,6]);
        else
            data_stocca = data_stocca([4,5,3,2]);
            sem = sem([4,5,3,2]);
            gstruct.ugroupdef = {gstruct.ugroupdef{1}(2:5)};
            gstruct.col = gstruct.col(2:5);
        end
        
        figure;%(fig(5))
        [~,~,barwidth] = barme([],data_stocca,sem,gstruct);
        ylabel('Ca^2^+ decay constant [ms]')
        ylim([0 300])
        FontResizer
        FigureResizer(5,8,[barwidth,0.4])
        tprint(fullfile(targetfolder_results,t2n_expcat('Fig.4-CaDecay',nneuron{n}.experiment)),'-pdf');
        
        data_stocca= [194,mean(nanmean(maxcai,3),1)*1e6,mean([0.91,1.16]*1000)];
        sem = [23,std(nanmean(maxcai,3),[],1)/sqrt(numel(tree))*1e6,std([0.91,1.16]*1000)/sqrt(2)];
        if ostruct.usemorph >= 4 && ~ostruct.newborn  % rat mature..add data
            data_stocca = data_stocca([1,4,5,3,2,6]);
            sem = sem([1,4,5,3,2,6]);
        else
            data_stocca = data_stocca([4,5,3,2]);
            sem = sem([4,5,3,2]);
        end
        
        figure;%(fig(6))
        [~,~,barwidth] = barme([],data_stocca,sem,gstruct);
        ylabel('Peak Ca amplitude [nM]')
        ax = gca;
        
        FontResizer
        if ax.YLim(2) > 2000
            ylim([0 3.5E+5])
        else
            ylim([0 1200])
        end
        FigureResizer(5,8,[barwidth,0.4])
        tprint(fullfile(targetfolder_results,t2n_expcat('Fig.4-CaAmp',nneuron{n}.experiment)),'-pdf');
        
    end
    if ~all(spiked)
        display('CAUTION: Not all cells spiked!')
    end
    fprintf('Mean Calcium decay time: Axon: %g +- %g nM Soma: %g +- %g ms Proximal: %g +- %g ms Distal: %g +- %g ms',[mean(nanmean(tw,3),1);std(nanmean(tw,3),[],1)])
    fprintf('Mean Calcium peak amplitude was: Axon: %g +- %g nM Soma: %g +- %g nM Proximal: %g +- %g nM Distal: %g +- %g nM',[mean(nanmean(maxcai,3),1);std(nanmean(maxcai,3),[],1)]*1e6)
    fprintf('Mean voltage attenuation @ %d µm: %g +- %g %% (s.e.m.)\n',thisdist,mean(bAPrelthisdist)*100,std(bAPrelthisdist)/sqrt(numel(tree))*100)
    if ostruct.show && ~isnan(data(1))
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