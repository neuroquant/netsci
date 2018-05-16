cd('+Cliquer');

strDir = fileparts(which('FindAll.m')); %[fileparts(tmp) '/cliquer/'];
% clean old stuff
system(['make -C ' 'cliquer/ clean']);
system(['make -C ' 'cliquer/']);

mex -v -Icliquer FindAll.c  cliquer/cliquer.o cliquer/graph.o cliquer/reorder.o

cd ..
