function fig = t2n_FIplot(targetfolder_data,neuron,ostruct,targetfolder_results)
% This function checks the neuron structure for correct definition of the
% used morphologies and returns info about it

% INPUTS
% targetfolder_data     folder which was given to t2n_currSteps, where the 
%                       data of the simulation lies
% neuron                t2n neuron structure (see documentation)
% ostruct               structure with fields defining some output
%                           figurewidth     width of figure to be created
%                           figureheigth    height of figure to be created
%                           savename        prefix filename of figures when saved
%                           ampprop         amplitude for which an extra
%                                           figure will be made and the maximal dV is
%                                           calculated
% targetfolder_results  folder where pdfs from figures should be saved. If
%                       not provided, figures will only be plotted
%
% OUTPUTS
% fig           figure handles to the output figures
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016, 2017 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if ~exist('targetfolder_results','var')
    targetfolder_results = [];
end
if ~isfield(ostruct,'duration') && ostruct.dataset < 7
    ostruct.duration = 200;
end

load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'),'voltVec','timeVec','numspikes','cstepsSpikingModel','tree','nneuron');

if isfield(ostruct,'color')
    modelcol = ostruct.color;
else
    modelcol = [1 0 0];
end

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(2))
    fig(1) = ostruct.handles(1);
    figure(fig(1))
else
    fig(1) = figure;hold all,
end
errorbar(cstepsSpikingModel*1000,mean(numspikes,1),std(numspikes,[],1),'Color',modelcol)
xlabel('Current step [pA]')
ylabel('Number of spikes')
FontResizer
if isfield(ostruct,'figureheight') && isfield(ostruct,'figurewidth')
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)
end
if ~isempty(targetfolder_results)
    if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
    else
        tprint(t2n_catName(targetfolder_results,'FI',neuron.experiment),'-pdf');
    end
end

if isfield(ostruct,'handles') && ~isempty(ostruct.handles) && ishandle(ostruct.handles(2))
    fig(2) = ostruct.handles(2);
    figure(fig(2))
else
    fig(2) = figure;hold all,
end
errorbar(cstepsSpikingModel*1000,mean(numspikes,1)/ostruct.duration * 1000,std(numspikes,[],1)/ostruct.duration * 1000,'Color',modelcol)
xlabel('current step [pA]')
ylabel('frequency [Hz]')
FontResizer
if isfield(ostruct,'figureheight') && isfield(ostruct,'figurewidth')
    FigureResizer(ostruct.figureheight,ostruct.figurewidth)
end
if ~isempty(targetfolder_results)
    if isfield(ostruct,'savename')  && ~isempty(ostruct.savename)
        if ~isempty(ostruct.savename)
            tprint(fullfile(targetfolder_results,[ostruct.savename,'_Hz']),'-pdf');
        end
    else
        tprint(t2n_catName(targetfolder_results,'FI_Hz',neuron.experiment),'-pdf');
    end
end