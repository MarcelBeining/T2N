function fig = t2n_plotCurrSteps(targetfolder_data,targetfolder_results,neuron,params,ostruct,steps)

if nargin < 6
    steps = [0.03,0.075]; % 30 and 75 pA
end

if nargin < 5
    ostruct.dataset = 2;
end
if ~isfield(ostruct,'show')
    if ostruct.usemorph >= 4  % rat
        ostruct.show = 2;
    else
        ostruct.show = 1:2;
    end
end
if numel(ostruct.amp)== 1
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
end

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;hold all,
end

if any(ostruct.show == 1)
    for f = 1:size(exp_iclamp,2)
        if numel(ostruct.amp)==1
            s = find(cstepsSpiking == ostruct.amp/1000);
            if ~isempty(s)
                hold all
                %         if params.realv
                this = exp_iclamp(:,f,s) - params.LJP ;   % corrected!
                if ostruct.dataset < 6
                    ylim([-70 -30]-round(params.LJP))
                end
                ylabel('Cell Voltage [mV] (corrected!)')
                title(sprintf('IClamp % 4.4g pA',cstepsSpiking(s)*1000));
                if exist('shoulder','var') && shoulder(f)
                    plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(this),'r');%exp_iclamp_mature(:,f,s)))
                else
                    plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(this),'Color',[0.4 0.4 0.4]);%exp_iclamp_mature(:,f,s)))
                end
                xlim([0 350])
            end
        else
            for s = 1:size(exp_iclamp,3)
                subplot(round(sqrt(size(exp_iclamp,3))),ceil(sqrt(size(exp_iclamp,3))),s)
                hold all
                if ostruct.dataset < 6
                    ylim([-70 -30]-round(params.LJP))
                end
                ylabel('Cell Voltage [mV] (corrected!)')
                title(sprintf('IClamp % 4.4g pA',cstepsSpiking(s)*1000));
                if ~all(isnan(exp_iclamp(:,f,s)))
                    %         if params.realv
                    this = exp_iclamp(:,f,s) - params.LJP ;   % corrected!
                    
                    %         else
                    %             this = exp_iclamp_mature(:,f,s) ;   % uncorrected!
                    %             ylim([-70 -30])
                    %             ylabel('Cell Voltage [mV] (uncorrected!)')
                    %             title(sprintf('IClamp % 4.4g pA',cstepsSpiking(s)*1000));
                    %         end
                    if exist('shoulder','var') && shoulder(f)
                        plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(this),'r');%exp_iclamp_mature(:,f,s)))
                    else
                        plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(this),'Color',[0.4 0.4 0.4]);%exp_iclamp_mature(:,f,s)))
                    end
                    
                    xlim([0 350])
                end
            end
        end
    end
    
    if ~all(isnan(exp_iclamp))
        Rin = (max(exp_iclamp(:,:,cstepsSpiking==0.01),[],1)-mean(exp_iclamp(1:50*rate,:,cstepsSpiking==0.01),1))/0.01;
        fprintf('\nMean Rin in Exp(@+10pA) is %g +- %g MOhm (s.e.m., +10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
    end
end

Rin = (cellfun(@(x) max(x),voltVec(:,cstepsSpikingModel==0.01))-cellfun(@(x) x(1),voltVec(:,cstepsSpikingModel==0.01)))/0.01;
fprintf('Mean Rin in Model(@+10pA) is %g +- %g MOhm (s.e.m., -10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))

if any(ostruct.show == 1)
    [~,ia,ib] = intersect(round(cstepsSpiking*1000),round(cstepsSpikingModel*1000)); % to account for if the steps were not the same in both
    nnan = ~isnan(cstepsSpikingModel);
    % if ~isempty(ia) && numel(cstepsSpikingModel) <= numel(cstepsSpiking)
    %     ia(nnan) = ia;
    %     ia(~nnan) = NaN;
    % end
    if ~isempty(ib)
        cstepsSpikingModel = cstepsSpikingModel(ib);
        voltVec = voltVec(:,ib);
        timeVec = timeVec(:,ib);
    end
else
    nnan = ones(numel(cstepsSpikingModel),1);
end
% if size(voltVec,2) > 25
%     voltVec = voltVec(:,1:25);
% end
if any(ostruct.show == 2)
    for f=1:size(voltVec,1)
        
        for s = 1:size(voltVec,2)
            if numel(ostruct.amp)==1
                hold all
                if ~isempty(timeVec{f,s})
                    plot(timeVec{f,s},squeeze(voltVec{f,s}),'LineWidth',1,'Color',tree{f}.col{1})
                end
            else
                if nnan(s)
                    if any(ostruct.show == 1)
                        subplot(round(sqrt(size(exp_iclamp,3))),ceil(sqrt(size(exp_iclamp,3))),ia(s))
                    else
                        subplot(round(sqrt(size(voltVec,2))),ceil(sqrt(size(voltVec,2))),s)
                    end
                    hold all
                    if ~isempty(timeVec{f,s})
                        plot(timeVec{f,s},squeeze(voltVec{f,s}),'LineWidth',1,'Color',tree{f}.col{1})
                    end
                end
            end
            %         ylim([-70 -30]-params.LJP)
            %         xlim([0 350])
        end
    end
end
if any(ostruct.show == 1) && ~isempty(intersect(cstepsSpiking,steps))
    if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && numel(ostruct.handles) > 1 && ishandle(ostruct.handles(2))
        fig(2) = ostruct.handles(2);
        figure(fig(2))
    else
        fig(2) = figure;hold all,
    end
    
    for ss = 1:2
        s = find(cstepsSpiking == steps(ss));
        if ~isempty(s)
            %             ff = [1 1 1];
            %             while numel(ff)~=numel(unique(ff))
            %                 ff=ceil(rand(3,1)*size(exp_iclamp,2));%1:size(exp_iclamp_mature,2)
            %             end
            
            if ostruct.newborn && ostruct.dataset == 2.28;
                ff = [2,7,12];
            else
                ff = [4 1 8];
            end
            for f=1:3%1:size(voltVec,1)
                %         if params.realv
                this = exp_iclamp(:,ff(f),s) - params.LJP ;   % corrected!
                p(f)= plot(1/rate:1/rate:size(exp_iclamp,1)/rate,squeeze(this),'LineWidth',1);%exp_iclamp_mature(:,f,s)))
            end
                set(p(1),'Color','k')
                set(p(2),'Color',[0.3 0.3 0.3])
                set(p(3),'Color',[0.6 0.6 0.6])
            
        end
    end
end
if any(ostruct.show == 2) && ~isempty(intersect(cstepsSpikingModel,steps))
    if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && numel(ostruct.handles) > 2  && ishandle(ostruct.handles(3))
        fig(3) = ostruct.handles(3);
        figure(fig(3))
    else
        fig(3) = figure;hold all,
    end
    for s = 1:numel(steps)
        if any(cstepsSpikingModel==steps(s))
            switch ostruct.usemorph
                case 1
                    if ostruct.newborn
                        ff = [3 4 7];
                    else
                        ff = 1:3;
                    end
                    if ostruct.reducecells
                        ff = 1;
                    end
                case 2
                    ff = [9 11 15];
                otherwise
                    ff = 1:3;
            end
            for f = ff%numel(tree)
                plot(timeVec{f,cstepsSpikingModel==steps(s)},squeeze(voltVec{f,cstepsSpikingModel==steps(s)}),'LineWidth',1,'Color',tree{f}.col{1})
            end
%             'g'
        end
    end
end
if any(ostruct.show == 3)
    if any(ostruct.show == 1)
        figure(fig(2))
        FontResizer
        FigureResizer(ostruct.figureheight,ostruct.figurewidth)
        xlabel('Time [ms]')
        ylabel('Membrane voltage [mV]')
        xlim([0 350])
        ylim([-85 -30])
        tprint(fullfile(targetfolder_results,'Fig.2-SpikingExp'),'-HR-pdf');
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
            tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-SpikingExp')),'-pdf');
        else
            tprint(fullfile(targetfolder_results,strcat('Fig.2-SpikingExp_',neuron.experiment)),'-pdf');
        end
    end
    if any(ostruct.show == 2) && ~isempty(intersect(cstepsSpikingModel,steps))
        figure(fig(3))
        FontResizer
        FigureResizer(ostruct.figureheight,ostruct.figurewidth)
        xlabel('Time [ms]')
        ylabel('Membrane voltage [mV]')
        if params.tstop > 350
            xlim([0 params.tstop])
        else
            xlim([0 350])
        end
        ylim([-85 -30])
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
            tprint(fullfile(targetfolder_results,strcat(ostruct.savename,'-SpikingModel')),'-pdf');
        else
            tprint(fullfile(targetfolder_results,strcat('Fig.2-SpikingModel_',neuron.experiment)),'-pdf');
        end
    end
end