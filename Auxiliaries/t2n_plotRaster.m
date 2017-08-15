function [] = t2n_plotRaster(varargin)
% This function creates a raster plot of a spike train. 
%
% INPUTS
% spikeMat      logical matrix with Ones where a spike occurs
%               alternatively it can be a vector with spiking times, or for
%               multiple spike trains, a cell array with each cell comprising
%               the spike time vectors
% tVec          corresponding time vector [ms]
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


switch nargin
    case 1
        % not used right now
    case 2
        spikeMat = varargin{1};
        tVec = varargin{2};
end
if ~islogical(spikeMat) % if spikeMat is not spikeMat but times of spike
    tmp = spikeMat;
    if ~iscell(tmp)
        tmp = {tmp};
    end
    spikeMat = false(numel(tmp),numel(tVec)); % init spikeMat
    for n = 1:numel(tmp)
        if iscell(tmp{n})
            tmp{n} = tmp{n}{1};
        end
        spikeMat(n,interp1(tVec,1:numel(tVec),tmp{n})) = true; %make spikeMat 1 at times
    end
end

hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 tVec(end)])
set(gca,'YTick',(1:trialCount))
set(gca,'YTickLabel',1:1:trialCount)
ylabel('Cell number')
xlabel('Time [ms]')