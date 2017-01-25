if(exist('lib/matlab-library','dir'))
	addpath(genpath('lib/matlab-library/lib/'));
	addpath(genpath('lib/matlab-library/test/'));
else
	error('matlab-library not installed at lib/')
end

	