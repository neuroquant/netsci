if(exist('lib/matlab-library','dir'))
	addpath('lib/matlab-library/lib/');
	addpath('lib/matlab-library/test/');
else
	error('matlab-library not installed at lib/')
end

addpath('netsci');
addpath(fullfile('netsci','utils')); 

import +community.*
import +tda.*

addpath('external/matlab-cliquer/')
addpath('external/k-clique/')
addpath('~/MATLAB/packages/matlab-bgl');
