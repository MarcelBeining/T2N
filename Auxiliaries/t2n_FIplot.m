function fig = t2n_FIplot(targetfolder_data,targetfolder_results,neuron,params,ostruct)
%
if ~isfield(ostruct,'duration') && ostruct.dataset < 7
    ostruct.duration = 200;
end
if ~isfield(ostruct,'show')
    ostruct.show = 1:2;
end

if numel(ostruct.amp)==1
    str = sprintf('_%dpA',ostruct.amp);
else
    str = '';
end
if params.cvode == 0
    load(expcat(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,str,'_fixed-dt')),'voltVec','timeVec','numspikes','params','cstepsSpikingModel','tree','nneuron')
else
    load(expcat(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,str)),'voltVec','timeVec','numspikes','params','cstepsSpikingModel','tree','nneuron')
end
if any(ostruct.usemorph == [2,3,5,6])  % artificial cells
    modelcol = [0 1 0];
else
    modelcol = [0 0 1];
end
[exp_iclamp,cstepsSpiking,rate] = load_ephys(ostruct.dataset,'CClamp');

if numel(ostruct.amp)==1
    s2= 1;
    s1 = find(cstepsSpiking == ostruct.amp/1000);
else
    if isfield(ostruct,'ampprop')
        s2=find(cstepsSpikingModel == ostruct.ampprop/1000);
        s1=find(cstepsSpiking == ostruct.ampprop/1000);
    else
        s2=find(cstepsSpikingModel == 0.09);
        s1 = find(cstepsSpiking == 0.09);
    end
end
if isempty(exp_iclamp)
    s1 = [];
end
instFImodel = NaN(numel(tree),floor(ostruct.duration/100));

if ~isempty(s2)
    for t = 1:numel(tree)
        for w = 1:floor(ostruct.duration/100)
            instFImodel(t,w) = sum(diff(voltVec{t,s2} > 0,1,1) == -1 & timeVec{t,s2}(2:end)>=55+(w-1)*100 & timeVec{t,s2}(2:end)<55+(w)*100)/0.1;
        end
    end
    if ~isempty(s1)
        tvec = (1/rate:1/rate:size(exp_iclamp,1)/rate)';
        for t = 1:size(exp_iclamp,2)
            for w = 1:floor(ostruct.duration/100)
                instFIexp(t,w) = sum(diff(exp_iclamp(:,t,s1) > 0,1,1) == -1 & tvec(2:end)>=55+(w-1)*100 & tvec(2:end)<55+(w)*100)/0.1;
            end
        end
    end
    
    if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
        fig(1) = ostruct.handles(1);
        figure(fig(1))
    else
        fig(1) = figure;hold all,
    end
    xvec = (1:floor(ostruct.duration/100)) * 100;
    if ~isempty(s1)
        hp = patch ([xvec (fliplr (xvec))], [(mean(instFIexp,1) + std(instFIexp,[],1)) (fliplr (mean(instFIexp,1) - std(instFIexp,[],1)))], [0 0 0]);
        hold on
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (xvec, mean(instFIexp,1), 'k')
    end
    errorbar(xvec,mean(instFImodel,1),std(instFImodel,[],1),'Color',modelcol)
    xlim([0 xvec(end)+100])
    set(gca,'XTick',[0,xvec])
    yl = get(gca,'YLim');
    ylim([0 yl(2)])
    xlabel('Intervals [ms]')
    ylabel('Instantaneous frequency [Hz]')
end




if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(2))
    fig(2) = ostruct.handles(2);
    figure(fig(2))
else
    fig(2) = figure;hold all,
end
if any(ostruct.show == 1)
    if ostruct.duration == 200 && ostruct.dataset ~= 0
        x = cstepsSpiking*1000;
        mFI = mean (squeeze((sum(diff(exp_iclamp > 0,1,1) == -1,1))));
        stdFI = std (squeeze((sum(diff(exp_iclamp > 0,1,1) == -1,1))));
        mFI = mFI(1:numel(x));
        stdFI= stdFI(1:numel(x));
        hp = patch ([x (fliplr (x))], [(mFI + stdFI) (fliplr (mFI - stdFI))],[0.7 0.7 0.7],'edgecolor','none' );
        if all(ostruct.show == 1)
        plot (x, mFI, 'k')
        end
    elseif ostruct.duration == 900 && exist(fullfile(params.path,'raw data','FI_Brenner.csv'),'file')
        dataBrenner = importdata(fullfile(params.path,'raw data','FI_Brenner.csv'));
        hp = patch ([50:50:400 (fliplr (50:50:400))], [(dataBrenner.data(2:2:end))' (fliplr (2*dataBrenner.data(1:2:end)'-dataBrenner.data(2:2:end)'))], [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (50:50:400, dataBrenner.data(1:2:end), 'k')
    elseif ostruct.duration == 1000 && exist(fullfile(params.path,'raw data','FI_Mehranfard15b_RAT.csv'),'file')
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard15b_RAT.csv'));
        MAstd = dataMA(7:end,2)-dataMA(1:6,2);
        hp = patch ([dataMA(1:6,1); (flipud (dataMA(1:6,1)))], [dataMA(1:6,2)-MAstd ;(flipud (dataMA(1:6,2)+MAstd))], [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (dataMA(1:6,1), dataMA(1:6,2), 'k')
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard14_RAT.csv'));
        MAstd = -dataMA.data(1:2:end-1)+dataMA.data(2:2:end);
        hp = patch ([(100:50:250)'; (fliplr (100:50:250)')], [dataMA.data(1:2:end-1)-MAstd ;(flipud (dataMA.data(1:2:end-1)+MAstd))], [0 0 0]);
        set (hp, 'facecolor',[0.7,0.7,0.7],'facealpha', 0.2, 'edgecolor', 'none')
        plot (100:50:250, dataMA.data(1:2:end-1), 'Color',[0.7,0.7,0.7])
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard15a_RAT.csv'));
        MAstd = dataMA(2:2:end,2)-dataMA(1:2:end-1,2);
        hp = patch ([(50:50:250)'; (flipud ((50:50:250)'))], [dataMA(1:2:end-1,2)-MAstd ;(flipud (dataMA(1:2:end-1,2)+MAstd))], [0 0 0]);
        set (hp, 'facecolor',[0.5,0.5,0.5],'facealpha', 0.2, 'edgecolor', 'none')
        plot (50:50:250, dataMA(1:2:end-1,2),  'Color',[0.5,0.5,0.5])
    end
end
if any(ostruct.show == 2) 
    errorbar(cstepsSpikingModel*1000,mean(numspikes,1),std(numspikes,[],1),'Color',modelcol)
end
if max(cstepsSpikingModel*1000) <= 120
    if ostruct.newborn
        ylim([0 10])
        set(gca,'YTick',0:2:10)
    else
        ylim([0 8])
        set(gca,'YTick',0:2:8)
    end
        xlim([0 120])
else
    xlim([40 310])
    ylim([0 60])
    set(gca,'YTick',0:20:60)
end
if isnan(ostruct.vmodel) && ostruct.usemorph >= 4 % AH99
    xlim([0 310])
    ylim([0 120])
    set(gca,'YTick',0:20:120)
end
xlabel('Current step [pA]')
ylabel('Number of spikes')
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf')
else
    tprint(fullfile(targetfolder_results,expcat('Fig.4-FI',neuron.experiment)),'-pdf')
end



if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(3))
    fig(3) = ostruct.handles(3);
    figure(fig(3))
else
    fig(3) = figure;hold all,
end
if any(ostruct.show == 1)
    if ostruct.duration == 200
        x = cstepsSpiking*1000;
        mFI = mFI / 200 * 1000; % number of spikes divided by time of stimulation (200 ms) = Hz
        stdFI = stdFI / ostruct.duration * 1000;
        hp = patch ([x (fliplr (x))], [(mFI + stdFI) (fliplr (mFI - stdFI))], [0.7 0.7 0.7],'edgecolor','none');
        hold on
        if all(ostruct.show == 1)
            plot (x, mFI, 'k')
        end
    elseif ostruct.duration == 900
        hp = patch ([50:50:400 (fliplr (50:50:400))], [(dataBrenner.data(2:2:end))' (fliplr (2*dataBrenner.data(1:2:end)'-dataBrenner.data(2:2:end)'))]/0.9, [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (50:50:400, dataBrenner.data(1:2:end)/0.9, 'k')
    elseif ostruct.duration == 1000 && exist(fullfile(params.path,'raw data','FI_Mehranfard15b_RAT.csv'),'file')
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard15b_RAT.csv'));
        MAstd = dataMA(7:end,2)-dataMA(1:6,2);
        hp = patch ([dataMA(1:6,1); (flipud (dataMA(1:6,1)))], [dataMA(1:6,2)-MAstd ;(flipud (dataMA(1:6,2)+MAstd))], [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        plot (dataMA(1:6,1), dataMA(1:6,2), 'k')
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard14_RAT.csv'));
        MAstd = -dataMA.data(1:2:end-1)+dataMA.data(2:2:end);
        hp = patch ([(100:50:250)'; (fliplr (100:50:250)')], [dataMA.data(1:2:end-1)-MAstd ;(flipud (dataMA.data(1:2:end-1)+MAstd))], [0 0 0]);
        set (hp, 'facecolor',[0.7,0.7,0.7],'facealpha', 0.2, 'edgecolor', 'none')
        plot (100:50:250, dataMA.data(1:2:end-1), 'Color',[0.7,0.7,0.7])
        dataMA = importdata(fullfile(params.path,'raw data','FI_Mehranfard15a_RAT.csv'));
        MAstd = dataMA(2:2:end,2)-dataMA(1:2:end-1,2);
        hp = patch ([(50:50:250)'; (flipud ((50:50:250)'))], [dataMA(1:2:end-1,2)-MAstd ;(flipud (dataMA(1:2:end-1,2)+MAstd))], [0 0 0]);
        set (hp, 'facecolor',[0.5,0.5,0.5],'facealpha', 0.2, 'edgecolor', 'none')
        plot (50:50:250, dataMA(1:2:end-1,2),  'Color',[0.5,0.5,0.5])
    end
end
if any(ostruct.show == 2)
    errorbar(cstepsSpikingModel*1000,mean(numspikes,1)/ostruct.duration * 1000,std(numspikes,[],1)/ostruct.duration * 1000,'Color',modelcol)
end
if ostruct.usemorph < 4 % mouse
    ylim([0 40])
    xlim([0 120])
else
    ylim([0 60])
    xlim([40 310])
end
if isnan(ostruct.vmodel) && ostruct.usemorph >= 4 % AH99
    xlim([0 310])
    ylim([0 120])
    set(gca,'YTick',0:20:120)
end
xlabel('current step [pA]')
ylabel('frequency [Hz]')
FontResizer
        FigureResizer(ostruct.figureheight,ostruct.figurewidth) 
if isfield(ostruct,'savename')  && ~isempty(ostruct.savename)
    if ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,[ostruct.savename,'_Hz']),'-pdf')
    end
else
    tprint(fullfile(targetfolder_results,expcat('Fig.4-FI_Hz',neuron.experiment)),'-pdf')
end