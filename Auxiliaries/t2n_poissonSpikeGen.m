function [spikeMat, tVec] = t2n_poissonSpikeGen(freq, par, nTrials)
% dt, tstop and tvec in ms, freq in Hz

if nargin < 3
    nTrials = 1;
end

nBins = floor(par.tstop/par.dt);
spikeMat = rand(nTrials, nBins) < freq*par.dt/1000;
tVec = 0:par.dt:par.tstop-par.dt;