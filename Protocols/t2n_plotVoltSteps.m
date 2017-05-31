function fig = t2n_plotVoltSteps(loadingfile,ostruct)
if nargin < 3
    ostruct = [];
end

steps = -120-12.1;

if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end

load(loadingfile,'currVec','vstepsModel','tree')

fig(1) = figure; hold all
fig(2) = figure;hold all
thiscurr = currVec;
if ostruct.subtract_hv
    for f = 1:numel(tree)%size(curr,1)
        for s = 1:size(thiscurr,2)
            thiscurr{f,s}(2,:) = thiscurr{f,s}(2,:) - mean(thiscurr{f,s}(2,thiscurr{f,s}(1,:)>=0 & thiscurr{f,s}(1,:)<=104));
        end
    end
end

str =  '';
figure(fig(1))
for f = 1:numel(tree)
    for s = 1:size(thiscurr,2)
        subplot(floor(sqrt(size(thiscurr,2))),ceil(sqrt(size(thiscurr,2))),s)
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