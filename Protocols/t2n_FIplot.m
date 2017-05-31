function fig = t2n_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct)
%
if ~isfield(ostruct,'duration') && ostruct.dataset < 7
    ostruct.duration = 200;
end

load(t2n_expcat(targetfolder_data,'Exp_Spiking',neuron.experiment),'voltVec','timeVec','numspikes','params','cstepsSpikingModel','tree','nneuron');

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
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
    tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
else
    tprint(fullfile(targetfolder_results,t2n_expcat('FI',neuron.experiment)),'-pdf');
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
FigureResizer(ostruct.figureheight,ostruct.figurewidth)
if isfield(ostruct,'savename')  && ~isempty(ostruct.savename)
    if ~isempty(ostruct.savename)
        tprint(fullfile(targetfolder_results,[ostruct.savename,'_Hz']),'-pdf');
    end
else
    tprint(fullfile(targetfolder_results,t2n_expcat('FI_Hz',neuron.experiment)),'-pdf');
end