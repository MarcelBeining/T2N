function [maxdv,shoulder] = t2n_plotdV(targetfolder_data,targetfolder_results,neuron,params,ostruct)
checkthis = 90;
if numel(ostruct.amp)==1
    str = sprintf('_%dpA',ostruct.amp);
else
    str = '';
end
if params.cvode == 0
    load(expcat(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,str,'_fixed-dt')))
else
    load(expcat(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,str)))
end

if any(ostruct.show == 1)
    [exp_iclamp,cstepsSpiking,rate] = load_ephys(ostruct.dataset,'CClamp');
    f = ones(1,20)/20;
    filt_exp_iclamp = filtfilt(f,1,exp_iclamp);
    dfilt = diff(exp_iclamp,1,1);
end

if  any(ostruct.usemorph == [2,3,5,6])  % artificial cells
    modelcol = [0 1 0];
else
    modelcol = [0 0 1];
end

fig(1) = figure('units','normalized','outerposition',[0 0 1 1]);ax(1) = axes;
if any(ostruct.show == 1)
    fig(2) = figure;ax(2) = axes;hold all
    xlabel('Cell Voltage [mV]')
    % ylabel('Rate of voltage change [mV/ms]')
    fig(3) = figure;ax(3) = axes;hold all
    xlabel('Cell Voltage [mV]')
    % ylabel('Rate of voltage change [mV/ms]')
end
fig(4) = figure;ax(4) = axes;hold all
xlabel('Cell Voltage [mV]')
% ylabel('Rate of voltage change [mV/ms]')
fig(5) = figure;ax(5) = axes;hold all
xlabel('Cell Voltage [mV]')
% ylabel('Rate of voltage change [mV/ms]')

if any(ostruct.show == 1)
    
    p = zeros(size(exp_iclamp,2),1);
    col = colorme(size(exp_iclamp,2),'-grk');
    figure(fig(1))
    for s = 1:size(exp_iclamp,3)
        
        subplot(floor(sqrt(size(exp_iclamp,3))),ceil(sqrt(size(exp_iclamp,3))),s)
        hold all
        for f=1:size(exp_iclamp,2)
            
            %         if params.realv
            this = exp_iclamp(2:end,f,s) - params.LJP ;   % corrected!
            %             ylim([-70 -30]-round(params.LJP))
            xlabel('Cell Voltage [mV] (corrected!)')
            %         else
            %             this = exp_iclamp(2:end,f,s) ;   % corrected!
            % %             ylim([-70 -30])
            %             xlabel('Cell Voltage [mV] (uncorrected!)')
            %         end
            ind = find(dfilt(:,f,s)>7,1,'first') + 5*rate;  % find first spike
            maxdv{1}(f,s) = max(dfilt(:,f,s))*rate;
            if isempty(ind)
                ind =1;
                linstyl = '-';
            else
                linstyl = ':';
            end
            dt = (diff(1/rate:1/rate:size(exp_iclamp,1)/rate,1))';
            
            plot(this(ind:end),dfilt(ind:end,f,s)./dt(ind:end),'LineWidth',1.5,'Color',colorme(col{f},'brighter'),'LineStyle',linstyl);%exp_iclamp(:,f,s)))  % rest of spikes is black
            p(f) = plot(this(1:ind),dfilt(1:ind,f,s)./dt(1:ind),'LineWidth',1.5,'Color',col{f},'LineStyle','-');%exp_iclamp(:,f,s)))  % mark first spike red
            if s == size(exp_iclamp,3)
                ind2 = find(this(1:ind)>= -23,1,'first');
                if ~isempty(ind2)
                    shoulder(f) = (dfilt(ind2,f,s)./dt(ind2)) > 75;
                end
            end
            if cstepsSpiking(s) == checkthis/1000
                %             figure(fig(2)), hold all
                plot(ax(3),this(ind:end),dfilt(ind:end,f,s)./dt(ind:end),'LineWidth',1.5,'Color',col{f},'LineStyle',linstyl);%exp_iclamp(:,f,s)))  % rest of spikes is dashed
                pp(f) = plot(ax(2),this(1:ind),dfilt(1:ind,f,s)./dt(1:ind),'LineWidth',1.5,'Color',col{f},'LineStyle','-');%exp_iclamp(:,f,s)))  % first spike is straight line
                %             figure(fig(1)), hold all
            end
            ylabel('dV')
            xlim([-80 80])          
            set(gca,'XTick',-80:40:80)
            ylim([-200 800])
            set(gca,'YTick',-200:200:800)
        end
        
        uistack((p),'top')
    end
else
    shoulder = NaN;
end

p = zeros(size(voltVec,1),1);
figure(fig(1))
for ss = 1:size(voltVec,2)
    for f=1:size(voltVec,1)
        if ostruct.usemorph >=4  % rat
            subplot(floor(sqrt(size(voltVec,2))),ceil(sqrt(size(voltVec,2))),ss)
        else
            subplot(floor(sqrt(size(exp_iclamp,3))),ceil(sqrt(size(exp_iclamp,3))),s)
        end
        hold all
        thisv = squeeze(voltVec{f,ss});
        maxdv{2}(f,ss) = max(diff(thisv,1,1))/params.dt;
        ind = find(voltVec{f,ss}>7,1,'first');  % find first

        thist = squeeze(timeVec{f,ss});
        if isempty(ind)
            ind =1;
            linstyl = '-';
            
        else
            if isnan(ostruct.vmodel)  %AH99 model.. has broader spikes...
                ind = find(thist >= thist(ind)+8,1,'first');% spike is finished within ~8ms
            else
                ind = find(thist >= thist(ind)+5,1,'first'); % spike is finished within ~5ms
            end
            linstyl = ':';
        end
        plot(thisv(ind+1:end),diff(thisv(ind:end),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle',linstyl);%,'Color',colorme(tree{f}.col{1},'brighter'))
        if ind ~=1
            p(f) = plot(thisv(2:ind),diff(thisv(1:ind),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle','-');
        end
        if cstepsSpikingModel(ss) == checkthis/1000
            %             figure(fig(3)),hold all
            plot(ax(5),thisv(ind+1:end),diff(thisv(ind:end),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle',linstyl);%,'Color',colorme(tree{f}.col{1},'brighter'))
            if ind ~=1
                p(f) = plot(ax(4),thisv(2:ind),diff(thisv(1:ind),1,1)./params.dt,'LineWidth',1.5,'Color',tree{f}.col{1},'LineStyle','-');
            end
        end
        xlim([-80 80])             
        set(gca,'XTick',-80:40:80)
            ylim([-200 800])
            set(gca,'YTick',-200:200:800)
%         end
    end
    uistack(setdiff(p,0),'top')
end

figure(fig(1))
if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
    tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotAll',ostruct.savename)),'-pdf')
else
    tprint(fullfile(targetfolder_results,expcat('PhasePlotAll',nneuron{1}.experiment)),'-pdf')
end
xlim([-80 80])             
set(gca,'XTick',-80:40:80)
    ylim([-200 800])
    set(gca,'YTick',-200:200:800)


if any(ostruct.show == 1)
    figure(fig(2))  % exp first spike
    xlim([-80 80])             
    set(gca,'XTick',-80:40:80)
    ylim([-200 800])
    set(gca,'YTick',-200:200:800)
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
    if  isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotExp',ostruct.savename)),'-pdf')
    else
        tprint(fullfile(targetfolder_results,expcat('PhasePlotExp',nneuron{1}.experiment)),'-pdf')
    end
    
    figure(fig(3))  %EXP 2nd spike
    xlim([-80 80])             
    set(gca,'XTick',-80:40:80)
    ylim([-200 800])
    set(gca,'YTick',-200:200:800)
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
    if  isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotExp2',ostruct.savename)),'-pdf')
    else
        tprint(fullfile(targetfolder_results,expcat('PhasePlotExp2',nneuron{1}.experiment)),'-pdf')
    end
end

if any(ostruct.show == 2)
    
    figure(fig(4))  
    xlim([-80 80])             
    set(gca,'XTick',-80:40:80)

        ylim([-200 800])
        set(gca,'YTick',-200:200:800)
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
    ylabel('')
    if  isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotModel',ostruct.savename)),'-pdf')
    else
        tprint(fullfile(targetfolder_results,expcat('PhasePlotModel',nneuron{1}.experiment)),'-pdf')
    end
    figure(fig(5))
    xlim([-80 80])             
    set(gca,'XTick',-80:40:80)
        ylim([-200 800])
        set(gca,'YTick',-200:200:800)

    ylabel('')
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
    if  isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,sprintf('%s-PhasePlotModel2',ostruct.savename)),'-pdf')
    else
        tprint(fullfile(targetfolder_results,expcat('PhasePlotModel2',nneuron{1}.experiment)),'-pdf')
    end
end
figure;hold all
if any(ostruct.show == 1)
     patch ([cstepsSpiking*1000 (fliplr (cstepsSpiking*1000))], [(mean(maxdv{1},1) + std(maxdv{1},[],1)) (fliplr (mean(maxdv{1},1) - std(maxdv{1},[],1)))], [0.6 0.6 0.6],'edgecolor','none');
        hold on
        plot (cstepsSpiking*1000, mean(maxdv{1},1), 'k')
end
if any(ostruct.show == 2)
    errorbar(cstepsSpikingModel*1000,mean(maxdv{2},1),std(maxdv{2},[],1),'Color',modelcol);%/sqrt(size(maxdv{2},1)) )
end
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth,ostruct)
ylabel('Maximal dV/dt [mV/ms]')
xlabel('current steps [pA]')
% legend('Experiment','Model')
if  isfield(ostruct,'savename') && ~isempty(ostruct.savename)
    tprint(fullfile(targetfolder_results,sprintf('%s-MaxdV',ostruct.savename)),'-pdf')
else
    tprint(fullfile(targetfolder_results,expcat('MaxdV',nneuron{1}.experiment)),'-pdf')
end

fprintf('Max dv of experiment: %g +- %g mV/ms (s.e.m.)\n',mean(maxdv{1}(:,ostruct.ampprop/1000==cstepsSpiking)),std(maxdv{1}(:,ostruct.ampprop/1000==cstepsSpiking))/sqrt(size(maxdv{1},1)) )
fprintf('Max dv of model: %g +- %g mV/ms (s.e.m.)\n',mean(maxdv{2}(:,ostruct.ampprop/1000==cstepsSpikingModel )),std(maxdv{2}(:,ostruct.ampprop/1000==cstepsSpikingModel))/sqrt(size(maxdv{2},1)) )