function [Rin, tau, cap, Vrest] = t2n_passTests(neuron,tree,params,targetfolder,ostruct)
% caution! ostruct.recordnode is vector giving each tree the location where
% to record!

Rin = NaN(1,numel(tree));
tau = NaN(1,numel(tree));
cap = NaN(1,numel(tree));
Vrest = NaN(1,numel(tree));

if isfield(ostruct,'recordnode') && isnumeric(ostruct.recordnode) && numel(ostruct.recordnode) == numel(tree)
    recordnode = ostruct.recordnode;
else
    recordnode = ones(numel(tree),1);
end
if isfield(ostruct,'stimnode') && isnumeric(ostruct.stimnode) && numel(ostruct.stimnode) == numel(tree)
    stimnode = ostruct.stimnode;
else
    stimnode = ones(numel(tree),1);
end
if ~isfield(ostruct,'figureheight')
    ostruct.figureheight = [];
    ostruct.figurewidth = [];
end
if ~isfield(ostruct,'passtest')
    ostruct.passtest = 'Mongiat';
end
if isfield(ostruct,'show') && ~any(ostruct.show == 0)
    options = '-s';
else
    options = '';
end
params.prerun = 1500;
del = 100;

params.tstart = 0;

params.cvode=0;
switch ostruct.passtest
    case 'Mongiat'
        Vh = -70-12.1;  % LJP corrected
        dur = 100;
%         params.cvode = 1;
        amp = -10; %mV for cap and Rin Mongiat 2009
        for t = 1:numel(tree)
            neuron.pp{t}.SEClamp = struct('node',stimnode(t),'times',[-100, del, del+dur],'amp', [Vh Vh+amp Vh],'rs',15); 
            neuron.record{t}.SEClamp = struct('record','i','node',recordnode(t));
        end
    case 'Mongiat2'
        dur = 500;
        amp = -0.01  ;      % 10pA for tau...Mongiat 2009
        [hstep, Vrest] = t2n_findCurr(tree,params,neuron,-80-12.1);  % LJP corrected
        for t = 1:numel(tree)
            neuron.pp{t}.IClamp = struct('node',stimnode(t),'times',[-400,del,del+dur],'amp', [hstep(t), hstep(t)+amp hstep(t)]); %n,del,dur,amp
            neuron.record{t}.cell = struct('record','v','node',recordnode(t));
        end
    case 'Std'
        dur = 500;
        amp = -0.01  ;      % 10pA for tau...Mongiat 2009
        [hstep, Vrest] = t2n_findCurr(tree,params,neuron,-80);
        for t = 1:numel(tree)
            neuron.pp{t}.IClamp = struct('node',stimnode(t),'times',[-400,del,del+dur],'amp', [hstep(t), hstep(t)+amp hstep(t)]); %n,del,dur,amp
            neuron.record{t}.cell = struct('record','v','node',recordnode(t));
        end
    case 'Mehranfard'
        dur = 300;
        amp = -0.05  ;      % -50pA Meranfahrd 2015
        [hstep, Vrest] = t2n_findCurr(tree,params,neuron,-70);
        for t = 1:numel(tree)
            neuron.pp{t}.IClamp = struct('node',stimnode(t),'times',[-400,del,del+dur],'amp', [hstep(t), hstep(t)+amp hstep(t)]); %n,del,dur,amp
            neuron.record{t}.cell = struct('record','v','node',recordnode(t));
        end
    case 'SH'
        dur = 500;
        amp = -0.003  ;      % 3pA Schmidt Hieber 2007
        
        for t = 1:numel(tree)
            neuron.pp{t}.IClamp = struct('node',stimnode(t),'del',del,'dur', dur,'amp', amp); %n,del,dur,amp
            neuron.record{t}.cell = struct('record','v','node',recordnode(t));
        end
    case 'Brenner'
        dur = 1000;
        amp = -0.02  ;      % -20pA Brenner 2005 from holding pot of -80mV
        [hstep, Vrest] = t2n_findCurr(tree,params,neuron,-80);
        for t = 1:numel(tree)
            neuron.pp{t}.IClamp = struct('node',stimnode(t),'times',[-400,del,del+dur],'amp', [hstep(t), hstep(t)+amp hstep(t)]); %n,del,dur,amp
            neuron.record{t}.cell = struct('record','v','node',recordnode(t));
        end
end
params.tstop = 500+2*dur;
params.dt = 2;
[out, ~] = t2n(tree,params,neuron,'-q-d');
fitstart = del+dur+2;
if ~isempty(strfind(options,'-s'))
    figure, hold on
end
if amp > 0
    fu = @le;
else
    fu = @ge;
end

switch ostruct.passtest
    case 'Mongiat'
        for t = 1:numel(tree)
            I0 = mean(out.record{t}.SEClamp.i{recordnode(t)}(1:del/params.dt+1)); 
            is = out.record{t}.SEClamp.i{recordnode(t)}((del+dur)/params.dt+1);
            if ~isfield(ostruct,'capacitance') || ostruct.capacitance == 1
                y = out.record{t}.SEClamp.i{recordnode(t)}(sign(amp)*out.record{t}.SEClamp.i{recordnode(t)} > sign(amp)*is)-is;
                x = out.t(sign(amp)*out.record{t}.SEClamp.i{recordnode(t)} > sign(amp)*is);
                cap(t) = trapz(x,y)/amp * 1000;
                fprintf('\n\nCapacitance is %g pF\n\n',cap(t))
            end
            Rin(t) = amp/(is-I0); %MOhm mV/nA
            fprintf('\n\nRin: %g MOhm ,     tau:  ms\n\n',  Rin(t));%,tauexp(t));
            if ~isempty(strfind(options,'-s'))
                plot(out.t,out.record{t}.SEClamp.i{recordnode(t)},'Color',tree{t}.col{1},'LineWidth',1)
                xlim([0 del+dur+100])
            end
        end
                    

    case {'SH','Mongiat2','Brenner','Mehranfard','Std'}
        for t = 1:numel(tree)
            V0 = mean(out.record{t}.cell.v{recordnode(t)}(1:del/params.dt+1));  %mV  only works with prerun
            v0(t) = V0;
            
            Vs = out.record{t}.cell.v{recordnode(t)}((del+dur)/params.dt+1); %mV   sensitive to noise since no mean
            vs(t) = Vs;
            Rin(t) = (Vs-V0)/(amp);  % MOhm
            
            xend = find(fu(out.record{t}.cell.v{recordnode(t)}((del+dur)/params.dt+1:end) , (Vs-V0)*0.1+V0),1,'first');  % take trace from current release until decay has reached 10% of max amplitude
            [a,~] = polyfit(out.t((fitstart)/params.dt+1:(fitstart)/params.dt+xend),log(sign(amp)*out.record{t}.cell.v{recordnode(t)}((fitstart)/params.dt+1:(fitstart)/params.dt+xend)-sign(amp)*V0),1);
            if ~isempty(strfind(options,'-s'))
                yf = NaN(1,numel(out.t));
                yf((fitstart)/params.dt+1:end) = sign(amp)*exp(out.t((fitstart)/params.dt+1:end) * a(1) +a(2))+V0;
                hold all
                if numel(tree) == t
                    xlabel('Time [ms]')
                end
                plot(out.t,out.record{t}.cell.v{recordnode(t)},'Color',tree{t}.col{1},'LineWidth',1)
                plot(out.t,yf,'Color','k','LineWidth',2,'LineStyle','--')
                if ceil(numel(tree)/2) == t
                    ylabel('Membrane Potential [mV]')
                end
            end
           
            tau(t) = -1/a(1);

            fprintf('\n\nRin: %g MOhm ,     tau: %g ms\n\n',  Rin(t),tau(t));%,tauth(t));
        end
        if ~isempty(strfind(options,'-s'))
            linkaxes
            if amp < 0
                ylim([floor(min(vs)) ceil(max(v0))])
            else
                ylim([floor(min(v0)) ceil(max(vs))])
            end
        end
        
end
if ~isempty(strfind(options,'-s'))
    xlim([0 2000])
    FontResizer
    FigureResizer(ostruct.figureheight,ostruct.figurewidth,[],ostruct)
    tprint(fullfile(targetfolder,strcat(sprintf('PassMeasure_%s_',ostruct.passtest),neuron.experiment)),'-pdf');
end
