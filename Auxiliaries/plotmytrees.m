function h = plotmytrees(tree,targetfolder,col,ostruct)
% this function plots and saves each tree in a nice way, eg for using in Adobe
% Illustrator
% this function is part of the T2N package
% INPUT
% tree: TREES toolbox tree cell array
% targetfolder: target folder for output files
% col (optional): cell array with rgb color values for each tree
% ostruct (optional): option structure with possible fields 'show' (1 = colored plotting
% but not good for putting in Illustrator; 2 = noncolored plotting),
% 'savename' (prefix name for image files).
% OUTPUT
% h: figure handles to the figures
%
% Copyright by Marcel Beining <marcel.beining@gmail.com>

if nargin < 4 || ~isfield(ostruct,'show')
    ostruct.show = 2;
end
if ~isfield(ostruct,'savename')
    ostruct.savename = 'Trees';
end
if nargin < 3 || isempty(col)
    col = colorme(numel(tree));
end

for t = 1:numel(tree)
    h=figure;hold all
    ptree = tran_tree(rot_tree(tran_tree(tree{t}),[],'-m3dY'),[0 300 0]);
    ptree.D(ptree.D<2) = 2;
    if ostruct.show == 2
        hp = plot_tree(ptree,col{t},[],[],[],'-b1');
        set(hp,'edgecolor','None')
    else
        if isfield(tree{t},'col')
            plot_tree(ptree,tree{t}.col{1},[],[],[],'-b1');
        else
            plot_tree(ptree,col{t},[],[],[],'-b1');
        end
    end
    ylim([100 750])
    xlim([-250,250])%*t+200])
    % xlim([200 1200])
    axis off
    ostruct.image = 1;
    if ostruct.show
        FigureResizer(10,10,0,ostruct)
        tprint(fullfile(targetfolder,sprintf('%s_%02g',ostruct.savename,t)),'-SHR-eps')
        %     tprint(fullfile(targetfolder,ostruct.savename),'-SHR-pdf')
        %     tprint(fullfile(targetfolder,ostruct.savename),'-SHR-svg')
%         tprint(fullfile2(targetfolder,ostruct.savename),'-SHR-png')
    end
end
