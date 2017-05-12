function fig = t2n_IVplot(loadingfile,targetfolder_results,ostruct)

if nargin < 3
    ostruct.dataset = 2;
end
if ~isfield(ostruct,'single')
    ostruct.single = 0;
end
if ~isfield(ostruct,'bablock')
    ostruct.bablock = 0;
end
if ~isfield(ostruct,'subtract_hv')
    ostruct.subtract_hv = 0;
end
if ~isfield(ostruct,'extract_kir')
    ostruct.extract_kir = 0;
end

load(loadingfile,'mholding_current','neuron','holding_voltage','steadyStateCurrVec','currVec','params','vstepsModel','tree','LJP')

if any(ostruct.show == 1) && ostruct.dataset ~= 0
    [exp_vclamp,vsteps,rate] = load_ephys(ostruct.dataset,'VClamp',ostruct.extract_kir);
    
    
    tvec = 1/rate:1/rate:size(exp_vclamp,1)/rate;
    for t = 1:size(exp_vclamp,2)
        d = exp_vclamp(:,t,vsteps==-80+10);
        is = mean(d(tvec>190&tvec<204));
        i0 = mean(d(tvec<104));
        Rin2(t) = (10)/(is-i0)*1000;  % Rin (in MOhm) mV/pA
        x = tvec(d < is & tvec' > 100);
        y = d(d < is & tvec' > 100)-is;
        y = y(x<=x(1)+50);  % voltage step length of original capacitance measurement was only 50 ms hence cut vector thereafter
        x = x(x<=x(1)+50);
        cap2(t) = trapz(x,y)/(-10); % calculate capacitance as integral (charge) divided by amplitude of voltage step
        
        d = exp_vclamp(:,t,vsteps==-80-10);
        is = mean(d(tvec>190&tvec<204));
        i0 = mean(d(tvec<104));
        Rin(t) = (-10)/(is-i0)*1000;  % Rin (in MOhm) mV/pA
        x = tvec(d < is & tvec' > 100);
        y = d(d < is & tvec' > 100)-is;
        y = y(x<=x(1)+50);  % voltage step length of original capacitance measurement was only 50 ms hence cut vector thereafter
        x = x(x<=x(1)+50);
        cap(t) = trapz(x,y)/(-10); % calculate capacitance as integral (charge) divided by amplitude of voltage step
        
    end
end
amp = [-10,10];
for a = 1:2
    for t = 1:numel(tree)
        ind = find(vstepsModel+params.LJP == -80+amp(a));
        I0 = mean(currVec{t,ind}(2,currVec{t,ind}(1,:)<104));
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
if any(ostruct.show == 1) && ostruct.dataset ~= 0 && (ostruct.dataset ~= 2.28)  % dont use that VClamp dataset, as it had been done after spiking experiment
    fprintf('Mean Rin in Exp(@-90mV-LJP) is %g +- %g MOhm (s.e.m., -10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
    fprintf('Mean Rin in Exp(@-70mV-LJP) is %g +- %g MOhm (s.e.m., +10mV)\n',mean(Rin2),std(Rin2)/sqrt(numel(Rin2)))
end
if  any(cellfun(@(x) ~any(x(1,:)>190&x(1,:)<204),currVec(:,vstepsModel+params.LJP==-80-10)))
    Rin = 1000*(-10)./cellfun(@(x) mean(x(2,(x(1,:)>175&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel+params.LJP==-80-10));
else
    Rin = 1000*(-10)./cellfun(@(x) mean(x(2,(x(1,:)>190&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel+params.LJP==-80-10));
end
fprintf('\nMean Rin in Model(@-90mV-LJP) is %g +- %g MOhm (s.e.m., -10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
Rin = 1000*(+10)./cellfun(@(x) mean(x(2,(x(1,:)>190&x(1,:)<204)))-mean(x(2,(x(1,:)<104))),currVec(:,vstepsModel+params.LJP==-80+10));
fprintf('Mean Rin in Model(@-70mV-LJP) is %g +- %g MOhm (s.e.m., +10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
if any(ostruct.show == 1) && ostruct.dataset ~= 0 && ostruct.dataset ~= 2.28  % dont use that VClamp dataset, as it had been done after spiking experiment
    fprintf('\nMean capacitance in Exp(@-90mV-LJP) is %g +- %g pF (s.e.m. -10mV)',mean(cap),std(cap)/sqrt(numel(cap)))
    fprintf('\nMean capacitance in Exp(@-70mV-LJP) is %g +- %g pF (s.e.m. +10mV)\n',mean(cap2),std(cap2)/sqrt(numel(cap2)))
end
fprintf('\nMean capacitance in Model(@-90mV-LJP) is %g +- %g pF (s.e.m. -10mV)',mean(capm(1,:)),std(capm(1,:))/sqrt(numel(capm(1,:))))
fprintf('\nMean capacitance in Model(@-70mV-LJP) is %g +- %g pF (s.e.m. +10mV)\n',mean(capm(2,:)),std(capm(2,:))/sqrt(numel(capm(2,:))))

if any(ostruct.show == 1) && ostruct.dataset ~= 0
    vstepsreal = vsteps - LJP;
    
    indhvexp = (vstepsreal == holding_voltage);
end
indhvmod = (vstepsModel == holding_voltage);


if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(1))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;clf;hold all,
end
if any(ostruct.show == 1) && ostruct.dataset ~= 0
    if exist('delind','var')
        meas_curr = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,setdiff(1:size(exp_vclamp,2),delind),:),1));
        basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,setdiff(1:size(exp_vclamp,2),delind),:),1));
    else
        meas_curr = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,:,:),1));
        basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,:,:),1));
    end
end
if ostruct.subtract_hv
    if any(ostruct.show == 1) && ostruct.dataset ~= 0
        meas_curr = meas_curr - basl;
    end
    steadyStateCurrVec = steadyStateCurrVec - repmat(mholding_current,size(steadyStateCurrVec,1),1);
end


Kirind_model = vstepsModel+LJP <= -110;
otherind_model = vstepsModel+LJP >= -70 & vstepsModel+LJP <= -50;

Restind_model = vstepsModel+LJP >= -100 & vstepsModel+LJP <= -70;


gKirModel = zeros(size(steadyStateCurrVec,2),1);
if any(ostruct.show == 1) && ostruct.dataset ~= 0
    gKirReal = zeros(size(meas_curr,1),1);
    Kirind_Mongiat = vsteps <= -110;
    otherind_Mongiat = vsteps >= -70 & vsteps <= -50;
    for t = 1:size(meas_curr,1)
        tmp = polyfit(vsteps(Kirind_Mongiat),meas_curr(t,Kirind_Mongiat),1) - polyfit(vsteps(otherind_Mongiat),meas_curr(t,otherind_Mongiat),1) ;
        gKirReal(t) = tmp(1);
    end
end
for t =1:size(steadyStateCurrVec,2)
    tmp = polyfit(vstepsModel(Kirind_model)+LJP,steadyStateCurrVec(Kirind_model,t)',1) - polyfit(vstepsModel(otherind_model)+LJP,steadyStateCurrVec(otherind_model,t)',1) ;
    gKirModel(t) = tmp(1);
end

tmp = polyfit(vstepsModel(Restind_model)+LJP,steadyStateCurrVec(Restind_model,t)',1);
fprintf('g at rest: %.3g nS\n',tmp(1))

if any(ostruct.show==1) && ostruct.dataset ~= 0

    mIV = mean(meas_curr,1);
    if ostruct.single
        plot(vstepsreal,meas_curr,'Color','k')
        for n = 1:size(meas_curr,1)
            vrest(n) = interp1(meas_curr(n,:),vstepsreal,0) ;
        end
        sprintf('Mean Vrest %g +- %g mV\n',mean(vrest),std(vrest))
    else
        
        stdIV = std (meas_curr,1);
        if ostruct.newborn
            if ostruct.dataset == 2.28 && exist(fullfile(params.path,'raw data','IV_29dpi_S3.csv'),'file')
                meas_curr = importdata(fullfile(params.path,'raw data','IV_29dpi_S3.csv'));
                mIV = meas_curr(1:19,2)';
                stdIV = meas_curr(20:end,2)' - mIV;
            else
                mIV = NaN(1,numel(vstepsreal));
                stdIV = NaN(1,numel(vstepsreal));
            end
        end
        hp = patch ([vstepsreal (fliplr (vstepsreal))], [(mIV + stdIV) (fliplr (mIV - stdIV))], [0 0 0]);
        set (hp, 'facealpha', 0.2, 'edgecolor', 'none')
        if all(ostruct.show == 1)
            plot (vstepsreal, mIV,'Color',[0 0 0],'LineWidth',1,'LineStyle','-')
        end
    end
end
line(vstepsModel,zeros(1,numel(vstepsModel)),'LineStyle','--','Color',[0.5 0.5 0.5])

if any(ostruct.show == 2)
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
    
    xlabel('Holding Voltage [mV] corrected')
end
if ostruct.usemorph <= 3 %mouse
    xlim([-140 -60])
else
    xlim([-125 -60])
end
if ostruct.newborn
    if ostruct.usemorph < 3 %mouse
        ylim([-150 50])
    else
        ylim([-200 100])
    end
else
    if ostruct.usemorph <= 3 %mouse
        ylim([-400 100])
    else
        ylim([-400 200])
    end
end
ylabel('Measured Current [pA]')



%% add paper data
vstepsreal = vstepsModel;
mIV = mean(steadyStateCurrVec,2);

Brenner05 = [-80,-20,600,10]; %HV, pA step, MOhm, stdevMOhm
Mongiat09 = [-70-12.1,-10,224,7]; %HV, mV step, MOhm, stdevMOhm  LJP corrected
SH07 = [-80,-3,308,26]; %HV, pA step, MOhm, stdevMOhm

Mongiat09y = [-70-12.1,-10,519,30]; %HV, mV step, MOhm, stdevMOhm  LJP corrected
Piatti11y21 = [-70-12.1,-10,636,94]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 21 dpi , LJP unknown
Piatti11y28 = [-70-12.1,-10,667,67]; %HV, mV step, MOhm, stdevMOhm  LJP corrected  %28 dpi , LJP unknown
Yang15y22a = [-70-12.1,-10,990,470]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 22 dpi Ascl , LJP unknown
Yang15y22b = [-70-12.1,-10,1040,640]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 22 dpi Glast, LJP unknown
Yang15y25a = [-70-12.1,-10,840,400]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 25 dpi Ascl , LJP unknown
Yang15y25b = [-70-12.1,-10,790,450]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 25 dpi Glast, LJP unknown
Yang15y28a = [-70-12.1,-10,560,200]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 28 dpi Ascl , LJP unknown
Yang15y28b = [-70-12.1,-10,660,250]; %HV, mV step, MOhm, stdevMOhm  LJP corrected   % 28 dpi Glast, LJP unknown


SH04 = [-80,-5,232,78]; % lila %HV, mV step, MOhm, stdevMOhm  RAT!!!
Staley92a = [-85,-15,107,NaN]; % rot %HV, mV step, MOhm, stdevMOhm  LJP corrected  % rat!
Staley92b = [-85,+15,228,14.2]; %rot %HV, mV step, MOhm, stdevMOhm  LJP corrected % rat!
MA14 = [-62,-50,230,15]; %orange %HV, pA step, MOhm, stdevMOhm  LJP corrected % rat!
Mehranfard = [-70,-50,295.6,11.5]; %grün

Brunner14yRAT21 = [-80,-10,378,20]; %HV, pA step, MOhm, stdevMOhm   % 21 dpi
Brunner14yRAT28 = [-80,-10,358,13]; %HV, pA step, MOhm, stdevMOhm   % 26 dpi

col = colorme({'red','yellow','dim green','violett','cyan','pink'});
if isfield(ostruct,'savename')
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)   % has to come before arrows to avoid distortions
end
if ~ostruct.newborn && ~ostruct.bablock
    if ostruct.usemorph <= 3  % mouse morphologies
        a=arrow([Brenner05(1),0+interp1(vstepsreal,mIV,Brenner05(1))],[Brenner05(1) + Brenner05(3) * Brenner05(2)/1000,Brenner05(2)+interp1(vstepsreal,mIV,Brenner05(1))]);
        set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([Mongiat09(1), interp1(vstepsreal,mIV,Mongiat09(1))],[Mongiat09(1) + Mongiat09(2), Mongiat09(2)/Mongiat09(3)*1000+interp1(vstepsreal,mIV,Mongiat09(1))]);
        set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([SH07(1),interp1(vstepsreal,mIV,SH07(1))],[SH07(1) + SH07(3) * SH07(2)/1000,SH07(2)+interp1(vstepsreal,mIV,SH07(1))]);
        set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
    else                            % rat model
        a=arrow([Staley92a(1),interp1(vstepsreal,mIV,Staley92a(1))],[Staley92a(1) + Staley92a(2),Staley92a(2)/Staley92a(3)*1000+interp1(vstepsreal,mIV,Staley92a(1))]);
        set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([Staley92b(1),interp1(vstepsreal,mIV,Staley92b(1))],[Staley92b(1) + Staley92b(2),Staley92b(2)/Staley92b(3)*1000+interp1(vstepsreal,mIV,Staley92b(1))]);
        set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([MA14(1),interp1(vstepsreal,mIV,MA14(1))],[MA14(1) + MA14(3) * MA14(2)/1000 , MA14(2)+interp1(vstepsreal,mIV,MA14(1))]);
        set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([Mehranfard(1),interp1(vstepsreal,mIV,Mehranfard(1))],[Mehranfard(1) + Mehranfard(3) * Mehranfard(2)/1000 , Mehranfard(2)+interp1(vstepsreal,mIV,Mehranfard(1))]);
        set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        a=arrow([SH04(1), interp1(vstepsreal,mIV,SH04(1))],[SH04(1) + SH04(2), SH04(2)/SH04(3)*1000+interp1(vstepsreal,mIV,SH04(1))]);
        set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
    end
elseif ostruct.newborn && ~ ostruct.bablock
    if ostruct.usemorph <= 3  % mouse morphologies
        a=arrow([Mongiat09y(1), interp1(vstepsreal,mIV,Mongiat09y(1))],[Mongiat09y(1) + Mongiat09y(2), Mongiat09y(2)/Mongiat09y(3)*1000+interp1(vstepsreal,mIV,Mongiat09y(1))]);
        set(a,'FaceColor',col{2},'EdgeColor',col{2},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        
        %         a=arrow([Piatti11y21(1), interp1(vstepsreal,mIV,Piatti11y21(1))],[Piatti11y21(1) + Piatti11y21(2), Piatti11y21(2)/Piatti11y21(3)*1000+interp1(vstepsreal,mIV,Piatti11y21(1))]);
        %         set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        
        a=arrow([Piatti11y28(1), interp1(vstepsreal,mIV,Piatti11y28(1))],[Piatti11y28(1) + Piatti11y28(2), Piatti11y28(2)/Piatti11y28(3)*1000+interp1(vstepsreal,mIV,Piatti11y28(1))]);
        set(a,'FaceColor',col{3},'EdgeColor',col{3},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        %
        %         a=arrow([Yang15y22a(1), interp1(vstepsreal,mIV,Yang15y22a(1))],[Yang15y22a(1) + Yang15y22a(2), Yang15y22a(2)/Yang15y22a(3)*1000+interp1(vstepsreal,mIV,Yang15y22a(1))]);
        %         set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        %
        %         a=arrow([Yang15y22b(1), interp1(vstepsreal,mIV,Yang15y22b(1))],[Yang15y22b(1) + Yang15y22b(2), Yang15y22b(2)/Yang15y22b(3)*1000+interp1(vstepsreal,mIV,Yang15y22b(1))]);
        %         set(a,'FaceColor',col{4},'EdgeColor',col{4},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        %
        %         a=arrow([Yang15y25a(1), interp1(vstepsreal,mIV,Yang15y25a(1))],[Yang15y25a(1) + Yang15y25a(2), Yang15y25a(2)/Yang15y25a(3)*1000+interp1(vstepsreal,mIV,Yang15y25a(1))]);
        %         set(a,'FaceColor',col{5},'EdgeColor',col{5},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        %
        %         a=arrow([Yang15y25b(1), interp1(vstepsreal,mIV,Yang15y25b(1))],[Yang15y25b(1) + Yang15y25b(2), Yang15y25b(2)/Yang15y25b(3)*1000+interp1(vstepsreal,mIV,Yang15y25b(1))]);
        %         set(a,'FaceColor',col{5},'EdgeColor',col{5},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        %
        a=arrow([Yang15y28a(1), interp1(vstepsreal,mIV,Yang15y28a(1))],[Yang15y28a(1) + Yang15y28a(2), Yang15y28a(2)/Yang15y28a(3)*1000+interp1(vstepsreal,mIV,Yang15y28a(1))]);
        set(a,'FaceColor',col{6},'EdgeColor',col{6},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        
        a=arrow([Yang15y28b(1), interp1(vstepsreal,mIV,Yang15y28b(1))],[Yang15y28b(1) + Yang15y28b(2), Yang15y28b(2)/Yang15y28b(3)*1000+interp1(vstepsreal,mIV,Yang15y28b(1))]);
        set(a,'FaceColor',col{6},'EdgeColor',col{6},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        
    else
        a=arrow([Brunner14yRAT21(1),0+interp1(vstepsreal,mIV,Brunner14yRAT21(1))],[Brunner14yRAT21(1) + Brunner14yRAT21(3) * Brunner14yRAT21(2)/1000,Brunner14yRAT21(2)+interp1(vstepsreal,mIV,Brunner14yRAT21(1))]);
        set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
        
        a=arrow([Brunner14yRAT28(1),0+interp1(vstepsreal,mIV,Brunner14yRAT28(1))],[Brunner14yRAT28(1) + Brunner14yRAT28(3) * Brunner14yRAT28(2)/1000,Brunner14yRAT28(2)+interp1(vstepsreal,mIV,Brunner14yRAT28(1))]);
        set(a,'FaceColor',col{1},'EdgeColor',col{1},'LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',0.5)
    end
end
%%
FontResizer

if isfield(ostruct,'savename')
    if ~isempty(ostruct.savename)
        
        tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
    end
else
    tprint(fullfile(targetfolder_results,expcat('Fig.2-IV',neuron.experiment)),'-pdf');
end
if any(ostruct.show == 1) && ostruct.dataset ~= 0
    fprintf('Kir Slope Conductance real cells %s\n',sprintf(' %.3g+-%.3g nS, ',mean(gKirReal),std(gKirReal)))
end
fprintf('Kir Slope Conductance model cells %s\n',sprintf(' %.3g+-%.3g nS, ',mean(gKirModel),std(gKirModel)))
