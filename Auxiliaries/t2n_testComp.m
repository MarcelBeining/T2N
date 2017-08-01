function tree = t2n_testComp

tree = struct('dA',sparse([0 0;1 0]),'X',[0;0],'Y',[0;0],'Z',[0;1],'name','testComp','D',[1;1],'R',[1;1],'rnames',{{'soma'}},'NID','testComp');
if ~exist(fullfile(pwd,'morphos'),'dir')
    warning('Caution! Folder "%s" does not exist in the current path and will be generated to save the hoc and mtr file of t2n_testComp',fullfile(pwd,'morphos'))
end
tree = t2n_writeTrees(tree,'',fullfile(pwd,'morphos','testComp.mtr'));