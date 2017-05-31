function [imp,freq] = t2n_resonance(amp,params,neuron,tree,ostruct,targetfolder_results)
if ~isfield(ostruct,'errorbar')
    ostruct.errorbar = 0;
end


tvec=0:params.dt:params.tstop;
freq = (0.5*(tvec+params.dt)/(1000));
vec = sin(2*pi*tvec.*freq/1000) * amp;
figure;plot(tvec/1000,freq)
xlabel('time [s]')
ylabel('frequency [Hz]')


figure;hold all;
xlabel('Frequency [Hz]')
ylabel('Impedance [M\Omega]')



hstep = find_curr(params,neuron,tree,ostruct.holding_voltage,[],'-q-d');
for t = 1:numel(tree)
    neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
    neuron.record{t}.cell = struct('node',1,'record','v');
end
% neuron_orig = neuron;
out = t2n(tree,params,neuron,'-q-d-w');
if any(cellfun(@(x) any(x.cell.v{1}>-30),out.record))
    warning('Caution! Spike was elicited!')
end


for t = 1:numel(tree)
    [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
    [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
    [locs,ia] = sort([locs;locs2]);
    pks = [pks;-pks2];
    pks = pks(ia)-ostruct.holding_voltage;
    imp{1}(:,t) = abs(pks/amp);
end
if ostruct.errorbar
    errorbar(freq(locs(2:2:end)-1),mean(imp{1}(2:2:end,:),2),std(imp{1}(2:2:end,:),[],2)/sqrt(numel(tree)),'color','k','linestyle','--');
    errorbar(freq(locs(1:2:end)-1),mean(imp{1}(1:2:end,:),2),std(imp{1}(1:2:end,:),[],2)/sqrt(numel(tree)),'color','k');
else
    plot(freq(locs(2:2:end)-1),mean(imp{1}(2:2:end,:),2),'color','k','linestyle','--');
end
% 
% % increase HCN
% for t = 1:numel(neuron.mech)
%     fields = fieldnames(neuron.mech{t});
%     for f1 = 1:numel(fields)
%         if isfield(neuron.mech{t}.(fields{f1}),'HCN')
%             neuron.mech{t}.(fields{f1}).HCN.gbar = neuron.mech{t}.(fields{f1}).HCN.gbar * 7;
%         end
% %         if isfield(neuron.mech{t}.(fields{f1}),'Kir21')
% %             neuron.mech{t}.(fields{f1}).Kir21.gkbar = neuron.mech{t}.(fields{f1}).Kir21.gkbar * 2;
% %         end
%     end
% end
% hstep = find_curr(params,neuron,tree,ostruct.holding_voltage,[],'-q-d');
% for t = 1:numel(tree)
%     neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
%     neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
% end
% out = t2n(tree,params,neuron,'-q-d-w');
% if any(out.record{t}.cell.v{1}>-30)
%     warndlg('Caution! Spike was elicited!')
% end
% for t = 1:numel(tree)
%     [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
%     [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
%     [locs,ia] = sort([locs;locs2]);
%     pks = [pks;-pks2];
%     pks = pks(ia)-ostruct.holding_voltage;
%     imp{2}(:,t) = abs(pks/amp);
% %     plot(freq(locs-1),imp,'b')
% end
% if ostruct.errorbar
%     errorbar(freq(locs(2:2:end)-1),mean(imp{2}(2:2:end,:),2),std(imp{2}(2:2:end,:),[],2)/sqrt(numel(tree)),'color','b','linestyle','--');
%     errorbar(freq(locs(1:2:end)-1),mean(imp{2}(1:2:end,:),2),std(imp{2}(1:2:end,:),[],2)/sqrt(numel(tree)),'color','b');
% else
%     plot(freq(locs(2:2:end)-1),mean(imp{2}(2:2:end,:),2),'color','r','linestyle','-');
% end
% 
% 
% % further increase Kir21
% for t = 1:numel(neuron.mech)
%     fields = fieldnames(neuron.mech{t});
%     for f1 = 1:numel(fields)
%         if isfield(neuron.mech{t}.(fields{f1}),'Kir21')
%             neuron.mech{t}.(fields{f1}).Kir21.gkbar = neuron.mech{t}.(fields{f1}).Kir21.gkbar * 2;
%         end
%     end
% end
% hstep = find_curr(params,neuron,tree,ostruct.holding_voltage,[],'-q-d');
% for t = 1:numel(tree)
%     neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
%     neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
% end
% out = t2n(tree,params,neuron,'-q-d-w');
% if any(out.record{t}.cell.v{1}>-30)
%     warndlg('Caution! Spike was elicited!')
% end
% for t = 1:numel(tree)
%     [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
%     [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
%     [locs,ia] = sort([locs;locs2]);
%     pks = [pks;-pks2];
%     pks = pks(ia)-ostruct.holding_voltage;
%     imp{3}(:,t) = abs(pks/amp);
% %     plot(freq(locs-1),imp,'r')
% end
% if ostruct.errorbar
%     errorbar(freq(locs(2:2:end)-1),mean(imp{3}(2:2:end,:),2),std(imp{3}(2:2:end,:),[],2)/sqrt(numel(tree)),'color','g','linestyle','--');
%     errorbar(freq(locs(1:2:end)-1),mean(imp{3}(1:2:end,:),2),std(imp{3}(1:2:end,:),[],2)/sqrt(numel(tree)),'color','g');
% else
%     plot(freq(locs(2:2:end)-1),mean(imp{3}(2:2:end,:),2),'color','b','linestyle','-');
% end
% 

neuron = t2n_blockchannel(neuron,'Kir21',100);
hstep = find_curr(params,neuron,tree,ostruct.holding_voltage,[],'-q-d');
for t = 1:numel(tree)
    neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
end
out = t2n(tree,params,neuron,'-q-d-w');
if any(cellfun(@(x) any(x.cell.v{1}>-30),out.record))
    warning('Caution! Spike was elicited!')
end
for t = 1:numel(tree)
    [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage,'NPeaks',450);
    [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage,'NPeaks',450);
    [locs,ia] = sort([locs;locs2]);
    pks = [pks;-pks2];
    pks = pks(ia)-ostruct.holding_voltage;
    imp{4}(:,t) = abs(pks/amp);
end
if ostruct.errorbar
    errorbar(freq(locs(2:2:end)-1),mean(imp{4}(2:2:end,:),2),std(imp{4}(2:2:end,:),[],2)/sqrt(numel(tree)),'color','g','linestyle','--');
    errorbar(freq(locs(1:2:end)-1),mean(imp{4}(1:2:end,:),2),std(imp{4}(1:2:end,:),[],2)/sqrt(numel(tree)),'color','g');
else
    plot(freq(locs(2:2:end)-1),mean(imp{4}(2:2:end,:),2),'color','b','linestyle','--');
end

neuron = t2n_blockchannel(neuron,'HCN',100);
hstep = find_curr(params,neuron,tree,ostruct.holding_voltage,[],'-q-d');
for t = 1:numel(tree)
    neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
end
out = t2n(tree,params,neuron,'-q-d-w');
if any(cellfun(@(x) any(x.cell.v{1}>-30),out.record))
    warning('Caution! Spike was elicited!')
end
for t = 1:numel(tree)
    [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
    [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',ostruct.holding_voltage);
    [locs,ia] = sort([locs;locs2]);
    pks = [pks;-pks2];
    pks = pks(ia)-ostruct.holding_voltage;
    imp{5}(:,t) = abs(pks/amp);
end
if ostruct.errorbar
    errorbar(freq(locs(2:2:end)-1),mean(imp{5}(2:2:end,:),2),std(imp{5}(2:2:end,:),[],2)/sqrt(numel(tree)),'color','g','linestyle','--');
    errorbar(freq(locs(1:2:end)-1),mean(imp{5}(1:2:end,:),2),std(imp{5}(1:2:end,:),[],2)/sqrt(numel(tree)),'color','g');
else
    plot(freq(locs(2:2:end)-1),mean(imp{5}(2:2:end,:),2),'color','r','linestyle','--');
end


FontResizer
FigureResizer(ostruct.figureheight,ostruct.figurewidth)

if ostruct.newborn
    str = '_newborn';
else
    str = '';
end
save(sprintf('%s%s%s_Resonance.mat',targetfolder_results,params.tname,str));

if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
    tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
else
    tprint(fullfile(targetfolder_results,sprintf('Fig5-Resonance_%s%s',params.tname,str)),'-pdf');
end
