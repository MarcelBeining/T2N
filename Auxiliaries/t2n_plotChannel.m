function t2n_plotChannel(neuron_orig,mcond_channel,region,outputFolder,hyper)
% This function activation and inactivation dynamics of an ion channel
%
% INPUTS
% neuron_orig       t2n neuron structure with already defined mechanisms
% mcond_channel     the name of the maximum conductance parameter of a
%                   channel, as it is in NEURON and NMODL (e.g. gbar_hh)
% region            name of the region in neuron_orig from which
%                   specifications for the channel should be taken from. If
%                   not provided, t2n takes the first specification it can
%                   find within a region
% outputFolder      (optional) folder where pdfs should be saved to
% hyper             boolean if channel is hyperpolarization activated
%                   (DEFAULT 0)
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************
if ~exist('hyper','var') || isempty(hyper)
    hyper = 0 ;
end
if ~exist('outputFolder','var') || isempty(outputFolder)
    outputFolder = '' ;
end

channel = strsplit(mcond_channel,'_'); mcond = channel{1};channel = channel{2};
cond_channel = strrep(mcond_channel,'bar','');

tree{1} = t2n_testComp;   % get the test compartment morphology

neuron.params = neuron_orig.params;
neuron.params.nseg = 1;
neuron.params.cvode = 0;
neuron.params.dt = 0.025;
neuron.mech{1}.all.(channel) = struct(mcond,0.1);

if exist('region','var') && ~isempty(region)  % get the ion channel specification from the region
    neuron.mech{1}.all.(channel) = neuron_orig.mech{1}.(region).(channel);
else   % else get it from the first region which incorporates the ion channel
    fields = fieldnames(neuron_orig.mech{1});
    for f = 1:numel(fields)
        if isfield(neuron_orig.mech{1}.(fields{f}),channel)
            neuron.mech{1}.all.(channel) = neuron_orig.mech{1}.(fields{f}).(channel);
        end
    end
end

neuron.record{1}.cell = struct('node',1,'record',cond_channel);

if hyper
    holding_voltage = 0;
    amp = -150:5:0;
else
    holding_voltage = -120;
    amp = -120:5:50;
end
[~,out] = t2n_VoltSteps(amp,[50 100 50],holding_voltage,neuron,tree);
maxG = zeros(numel(amp),1);
tauA = maxG;
figure;hold all;
xlabel('Time [ms]')
ylabel('Conductance of channel')
for a = 1:numel(amp)
    thisVec = out{a}.record{1}.cell.(cond_channel){1};
    [maxG(a),indmaxG] = max(thisVec);
    indStart = find(out{a}.t >= 50,1,'first');
    [ft,~] = polyfit(out{a}.t(indStart:indmaxG-1),log(-thisVec(indStart:indmaxG-1)+thisVec(indmaxG)),1);
    if imag(ft(1)) ~= 0
        ft = [NaN NaN];
    end
    plot(out{a}.t,thisVec)
    plot(out{a}.t(indStart:indmaxG),-exp(ft(1)*out{a}.t(indStart:indmaxG)+ft(2))+thisVec(indmaxG),'r--')
    
    tauA(a) = -1/ft(1);
end
maxG = maxG/max(maxG);
try
    ampHalf = interp1(maxG,amp,0.5);
catch
    ampHalf = NaN;
end
figure;plot(amp,maxG)
hold on;
scatter(ampHalf,0.5,'xr')
xlabel('Voltage [mV]')
ylabel('Norm. max. Conductance')
title('Activation curve')
xlim([amp(1) amp(end)])
ylim([0 1])
FontResizer
FigureResizer(5,5)
tprint(fullfile(outputFolder,sprintf('ActivCurve_%s_%+.3g',mcond_channel,ampHalf)),'-pdf-R')

tauA(maxG <= 0.0005) = NaN;    % delete nonsense fittings
figure;plot(amp,tauA)
xlabel('Voltage [mV]')
ylabel('activation time constant [ms]')
xlim([amp(1) amp(end)])
title('Activation time constant')
FontResizer
FigureResizer(5,5)
tprint(fullfile(outputFolder,sprintf('ActivTime_%s',mcond_channel)),'-pdf-R')
        

%% inactivation
% neuron.params.cvode = 1;
% neuron.params.dt = 0.025;

if hyper
    holding_voltage = [-150,-150];
    amp = -150:5:0;
else
    holding_voltage = [-120,50];
    amp = -120:5:50;
end
predur = 50;
dur = 1000;
actdur = 50;
[~,out] = t2n_VoltSteps(amp,[predur dur actdur],holding_voltage,neuron,tree);
maxG = zeros(numel(amp),1);
tauI = maxG;

figure;
hold all
xlabel('Time [ms]')
ylabel('Conductance of channel')
for a = 1:numel(amp)
    thisVec = out{a}.record{1}.cell.(cond_channel){1};
    maxG(a) = max(thisVec(out{a}.t>=predur+dur & out{a}.t<=predur+actdur+dur));
    [maxG1,indStart] = max(thisVec(out{a}.t>=predur & out{a}.t<=predur+actdur));
    indStart = indStart + find(out{a}.t>=predur,1,'first');
    indEnd = find(out{a}.t>=predur+dur,1,'first')-1;
    indEnd2 = find(thisVec(indStart:indEnd) < maxG1/25,1,'first') + indStart -1;
    if isempty(indEnd2)
        indEnd2 = indEnd;
    end
    [ft,~] = polyfit(out{a}.t(indStart:indEnd2-1),log(thisVec(indStart:indEnd2-1)-thisVec(indEnd2)),1);
    if imag(ft(1)) ~= 0
        ft = [NaN NaN];
    end
    plot(out{a}.t,thisVec)
    plot(out{a}.t(indStart:indEnd),exp(ft(1)*out{a}.t(indStart:indEnd)+ft(2))+thisVec(indEnd2),'r--')
    tauI(a) = -1/ft(1);
end
maxG = maxG/max(maxG);
try
    ampHalf = interp1(maxG,amp,0.5);
catch
    ampHalf = NaN;
end
figure;plot(amp,maxG)
hold on;
scatter(ampHalf,0.5,'xr')
xlabel('Voltage [mV]')
ylabel('Norm. max. Conductance')
title('Inactivation curve')
xlim([amp(1) amp(end)])
ylim([0 1])
FontResizer
FigureResizer(5,5)
tprint(fullfile(outputFolder,sprintf('InactivCurve_%s_%+.3g',mcond_channel,ampHalf)),'-pdf-R')

tauI(maxG >= 0.999) = NaN;   % delete nonsense fittings
figure;plot(amp,tauI)
xlabel('Voltage [mV]')
ylabel('Inactivation time constant [ms]')
xlim([amp(1) amp(end)])
title('Inactivation time constant')
FontResizer
FigureResizer(5,5)
tprint(fullfile(outputFolder,sprintf('InactivTime_%s',mcond_channel)),'-pdf-R')

