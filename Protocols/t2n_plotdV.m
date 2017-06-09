function [maxdv,fig] = t2n_plotdV(targetfolder_data,targetfolder_results,neuron,params,ostruct)
checkthis = 90;
load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'))

modelcol = [1 0 0];

fig(1) = figure('units','normalized','outerposition',[0 0 1 1]);ax(1) = axes;
fig(2) = figure;hold all
fig2(1) = figure;ax(4) = axes;hold all
xlabel('Cell Voltage [mV]')
fig2(2) = figure;ax(5) = axes;hold all
xlabel('Cell Voltage [mV]')

p = zeros(size(voltVec,1),1);
figure(fig(1))
for ss = 1:size(voltVec,2)
    for f=1:size(voltVec,1)
        subplot(floor(sqrt(size(voltVec,2))),ceil(sqrt(size(voltVec,2))),ss)
        hold all
        thisv = squeeze(voltVec{f,ss});
        maxdv{2}(f,ss) = max(diff(thisv,1,1))/params.dt;
        ind = find(voltVec{f,ss}>0,1,'first');  % find first AP
        
        thist = squeeze(timeVec{f,ss});
        if isempty(ind)
            ind =1;
            linstyl = '-';
            
        else
            ind = find(thist >= thist(ind)+5,1,'first'); % first spike should be finished within ~6ms
            linstyl = ':';
        end
        plot(thisv(ind+1:end),diff(thisv(ind:end),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle',linstyl);%,'Color',colorme(tree{f}.col{1},'brighter'))
        if ind ~=1
            p(f) = plot(thisv(2:ind),diff(thisv(1:ind),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle','-');
        end
        if cstepsSpikingModel(ss) == checkthis/1000
            plot(ax(5),thisv(ind+1:end),diff(thisv(ind:end),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle',linstyl);%,'Color',colorme(tree{f}.col{1},'brighter'))
            if ind ~=1
                p(f) = plot(ax(4),thisv(2:ind),diff(thisv(1:ind),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle','-');
            end
        end
        xlim([-80 80])
        set(gca,'XTick',-80:40:80)
        ylim([-200 800])
        set(gca,'YTick',-200:200:800)
    end
    uistack(setdiff(p,0),'top')
end

figure(fig(1))
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotAll',ostruct.savename)),'-pdf')


figure(fig2(1))
xlim([-80 80])
set(gca,'XTick',-80:40:80)

ylim([-200 800])
set(gca,'YTick',-200:200:800)
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
ylabel('')
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotModel',ostruct.savename)),'-pdf')
figure(fig2(2))
xlim([-80 80])
set(gca,'XTick',-80:40:80)
ylim([-200 800])
set(gca,'YTick',-200:200:800)

ylabel('')
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotModel2',ostruct.savename)),'-pdf')

figure(fig(2))
errorbar(cstepsSpikingModel*1000,mean(maxdv{2},1),std(maxdv{2},[],1),'Color',modelcol);%/sqrt(size(maxdv{2},1)) )
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
ylabel('Maximal dV/dt [mV/ms]')
xlabel('current steps [pA]')
tprint(fullfile(targetfolder_results,sprintf('%s-MaxdV',ostruct.savename)),'-pdf')


fprintf('Max dv of model: %g +- %g mV/ms (s.e.m.)\n',mean(maxdv{2}(:,ostruct.ampprop==cstepsSpikingModel )),std(maxdv{2}(:,ostruct.ampprop==cstepsSpikingModel))/sqrt(size(maxdv{2},1)) )