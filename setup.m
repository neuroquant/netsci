% Be sure to install github.com:TheEtkinLab/matlab-library in the lib/ directory. 

MATLABLIB='lib/matlab-library';

if(exist(MATLABLIB,'dir'))
	addpath(fullfile(MATLABLIB,'lib'));
	addpath(fullfile(MATLABLIB,'test'));
    addpath(fullfile(MATLABLIB,'lib','BrewerMap'));
else
	error('matlab-library not installed at lib/')
end

addpath('netsci');
addpath('data');
addpath(fullfile('netsci','utils')); 
addpath(fullfile('netsci','plot')); 

import +community.*
import +tda.*

addpath('external/matlab-cliquer/')
% addpath('external/k_clique/')
% addpath('~/MATLAB/packages/matlab-bgl');
