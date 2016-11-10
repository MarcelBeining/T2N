folder = uigetdir(pwd,'Please give a folder where the model structure should be initialized');
cd(folder)

mkdir(folder,'lib_mech')
mkdir(folder,'lib_genroutines')
mkdir(folder,'lib_custom')
mkdir(folder,'morphos')