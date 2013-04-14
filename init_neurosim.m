
% this initalization script sets up the execution search path

disp('Initializing neurosim environment ...');
clear all;


% assume baseDir is the location of this init script
global baseDir;
[baseDir, name, ext] = fileparts(mfilename('fullpath'));

% add include path for neurosim library
restoredefaultpath();                        % reset path to base
addpath( [ baseDir '/neurosim' ] );
addpath( [ baseDir '/util' ] );
addpath( [ baseDir '/sample_data' ] );

% set a flag so that we know initialization has been done
global neurosim_initialized; neurosim_initialized=true;

