function [spikeMat, tVec] = t2n_poissonSpikeGen(freq, params, nTrials)
% dt, tstop and tvec in ms, freq in Hz

if nargin < 3
    nTrials = 1;
end

%    params.dt = params.dt/1000; %ms to s 

nBins = floor(params.tstop/params.dt);
spikeMat = rand(nTrials, nBins) < freq*params.dt/1000;
tVec = 0:params.dt:params.tstop-params.dt;