function fig = t2n_IVplot(loadingfile,ostruct)

if ~isfield(ostruct,'single')
    ostruct.single = 0;
end
if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end

load(loadingfile,'mholding_current','neuron','holding_voltage','steadyStateCurrVec','currVec','vstepsModel','tree')

amp = [-10,10];
capm = zeros(2,numel(tree));
for a = 1:2
    for t = 1:numel(tree)
        ind = find(vstepsModel == holding_voltage+amp(a));
        %         I0 = mean(currVec{t,ind}(2,currVec{t,ind}(1,:)<104));
        if ~any(currVec{t,ind}(1,:)>190 & currVec{t,ind}(1,:)<204)
            is(t) = mean(currVec{t,ind}(2,currVec{t,ind}(1,:)>175 & currVec{t,ind}(1,:)<204));
        else
            is(t) = mean(currVec{t,ind}(2,currVec{t,ind}(1,:)>190 & currVec{t,ind}(1,:)<204));
        end
        y = currVec{t,ind}(2,sign(amp(a))*currVec{t,ind}(2,:) > sign(amp(a))*is(t))-is(t);
        x = currVec{t,ind}(1,sign(amp(a))*currVec{t,ind}(2,:) > sign(amp(a))*is(t));
        capm(a,t) = trapz(x,y)/amp(a);
    end
end
if  any(cellfun(@(x) ~any(x(1,:)>190&x(1,:)<204),currVec(:,vstepsModel==holding_voltage+amp(1))))
    Rin = 1000*(-10)./cellfun(@(x) mean(x(2,(x(1,:)>175&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel==holding_voltage+amp(1)));
else
    Rin = 1000*(-10)./cellfun(@(x) mean(x(2,(x(1,:)>190&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel==holding_voltage+amp(1)));
end
fprintf('\nMean Rin in Model(@%gmV) is %g +- %g MOhm (s.e.m., -10mV)\n',holding_voltage+amp(1),mean(Rin),std(Rin)/sqrt(numel(Rin)))
Rin = 1000*(+10)./cellfun(@(x) mean(x(2,(x(1,:)>190&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel==holding_voltage+amp(2)));
fprintf('Mean Rin in Model(@%gmV) is %g +- %g MOhm (s.e.m., +10mV)\n',holding_voltage+amp(2),mean(Rin),std(Rin)/sqrt(numel(Rin)))
fprintf('\nMean capacitance in Model(@%gmV) is %g +- %g pF (s.e.m. -10mV)',holding_voltage+amp(1),mean(capm(1,:)),std(capm(1,:))/sqrt(numel(capm(1,:))))
fprintf('\nMean capacitance in Model(@%gmV) is %g +- %g pF (s.e.m. +10mV)\n',holding_voltage+amp(2),mean(capm(2,:)),std(capm(2,:))/sqrt(numel(capm(2,:))))

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;clf;hold all,
end
if ostruct.subtract_hv
    steadyStateCurrVec = steadyStateCurrVec - repmat(mholding_current,size(steadyStateCurrVec,1),1);
end


Kirind_model = vstepsModel <= -121.1;
otherind_model = vstepsModel >= -82.1 & vstepsModel <= -62.1;
Restind_model = vstepsModel >= -112.1 & vstepsModel <= -82.1;
gKirModel = zeros(size(steadyStateCurrVec,2),1);
for t =1:size(steadyStateCurrVec,2)
    tmp = polyfit(vstepsModel(Kirind_model),steadyStateCurrVec(Kirind_model,t)',1) - polyfit(vstepsModel(otherind_model),steadyStateCurrVec(otherind_model,t)',1) ;
    gKirModel(t) = tmp(1);
end
tmp = polyfit(vstepsModel(Restind_model),steadyStateCurrVec(Restind_model,t)',1);
fprintf('Conductance at rest: %.3g nS\n',tmp(1))
fprintf('Kir Slope Conductance model cells %s\n',sprintf(' %.3g+-%.3g nS, ',mean(gKirModel),std(gKirModel)))

line(vstepsModel,zeros(1,numel(vstepsModel)),'LineStyle','--','Color',[0.5 0.5 0.5])

if ostruct.single
    p = plot(vstepsModel,steadyStateCurrVec);
    for t =1:size(steadyStateCurrVec,2)
        set(p(t),'color',tree{t}.col{1})
    end
elseif  any(ostruct.usemorph == [2,3,5,6])  % artificial cells
    errorbar(vstepsModel,mean(steadyStateCurrVec,2),std(steadyStateCurrVec,[],2)/sqrt(size(steadyStateCurrVec,2)),'Color',[0 1 0],'LineWidth',1);
else
    errorbar(vstepsModel,mean(steadyStateCurrVec,2),std(steadyStateCurrVec,[],2)/sqrt(size(steadyStateCurrVec,2)),'Color',[0 0 1],'LineWidth',1);
end

xlabel('Holding Voltage [mV]')
ylabel('Current [pA]')
