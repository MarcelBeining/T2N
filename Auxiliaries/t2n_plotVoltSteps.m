function t2n_plotVoltSteps(loadingfile,targetfolder_results,ostruct)
if nargin < 3
    ostruct = [];
end

steps = -120-12.1;

if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end
if ~isfield(ostruct,'show')
    ostruct.show = 1:2;
end

load(loadingfile,'mholding_current','neuron','holding_voltage','steadyStateCurrVec','currVec','params','vstepsModel','tree')
if any(ostruct.show==1)
    [exp_vclamp,vsteps,rate] = load_ephys(ostruct.dataset,'VClamp');
end

fig(1) = figure; hold all
fig(2) = figure;hold all
thiscurr = currVec;
if ostruct.subtract_hv
    if any(ostruct.show==1)
        basel = mean(exp_vclamp(0*rate+1:104*rate+1,:,:),1);
        exp_vclamp = exp_vclamp(:,:,:) - repmat(basel,size(exp_vclamp,1),1,1); % subtract current at baseline holding voltage (as Mongiat did)
    end
    for f = 1:numel(tree)%size(curr,1)
        for s = 1:size(thiscurr,2)
            thiscurr{f,s}(2,:) = thiscurr{f,s}(2,:) - mean(thiscurr{f,s}(2,thiscurr{f,s}(1,:)>=0 & thiscurr{f,s}(1,:)<=104));
        end
    end
end

str =  '';
figure(fig(1))
if any(ostruct.show == 1)
    for f = 1:size(exp_vclamp,2)
        for s = 1:size(exp_vclamp,3)
            subplot(floor(sqrt(size(exp_vclamp,3))),ceil(sqrt(size(exp_vclamp,3))),s)
            hold all
            plot(1/rate:1/rate:size(exp_vclamp,1)/rate,squeeze(exp_vclamp(:,f,s)))
            
            ylabel('Current [pA]')
            xlabel('Time [ms]')
            ylim([-200,200])
            title(sprintf('VClamp % 4.4g mV%s',vsteps(s),str));
            xlim([0 300])
            if any(vsteps(s) == steps)
                figure(fig(2))
                plot(1/rate:1/rate:size(exp_vclamp,1)/rate,squeeze(exp_vclamp(:,f,s)))
                ylabel('Current [pA]')
                xlabel('Time [ms]')
                ylim([-400,200])
                title(sprintf('VClamp % 4.4g mV%s',vsteps(s),str));
                xlim([0 300])
                figure(fig(1))
            end
        end
    end
end

if any(ostruct.show == 2)
    for f = 1:numel(tree)
        for s = 1:size(thiscurr,2)
            if any(ostruct.show == 1)
                subplot(floor(sqrt(size(exp_vclamp,3))),ceil(sqrt(size(exp_vclamp,3))),find(vstepsModel(s)==vsteps))
            else
                subplot(floor(sqrt(size(thiscurr,2))),ceil(sqrt(size(thiscurr,2))),s)
            end
            hold all
            plot(thiscurr{f,s}(1,:),thiscurr{f,s}(2,:),'LineWidth',3,'LineStyle','-','Color',tree{f}.col{1})
            ylabel('Current [pA]')
            xlabel('Time [ms]')
            ylim([-200,200])
            title(sprintf('VClamp % 4.4g mV%s',vstepsModel(s),str));
            if any(vstepsModel(s) == steps)
                figure(fig(2))
                plot(thiscurr{f,s}(1,:),thiscurr{f,s}(2,:),'LineWidth',1,'LineStyle','-','Color',tree{f}.col{1})
                ylabel('Current [pA]')
                xlabel('Time [ms]')
                ylim([-500,200])
                title(sprintf('VClamp % 4.4g mV%s',vstepsModel(s),str));
                xlim([0 300])
                figure(fig(1))
            end
        end
    end
end

figure(fig(2))
FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename')
    if ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
    end
else
    tprint(fullfile(targetfolder_results,expcat('Fig.2-IV_dyn',neuron.experiment)),'-pdf');
end