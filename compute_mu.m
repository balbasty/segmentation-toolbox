% TODO---------------------------------------------------------------------
% -Do not warp images beforehand
% -Push instead of interpolate
% -Change parameterisation of registration
% -Add regularisation for W and n update in update_pr (for more stability)
% -Speed up L calculation for different parts (e.g. use ll from mnom_slice in reg and bf)
% -Improve kmeans, maybe include regular GMM
% -Use only single again
% -Run on MRI+CT
% -Convolve (smooth) munum and muden
% -Regularise weight update (multiply 1 in the numerator)

% TEST---------------------------------------------------------------------
% -Try reg w/o weight updates

% Qs-----------------------------------------------------------------------
% -Should I get an estimate for the hidden responsibilities?

clear;

addpath('./core')
addpath('./util')

%--------------------------------------------------------------------------
% Initial parameters
%--------------------------------------------------------------------------

pars = [];

pars.imdir   = '/home/mbrud/Dropbox/PhD/Data/IXI-2D/';
pars.tempdir = './temp/'; % A temporary folder which will contain the preprocessed and warped images 

pars.N = 8; % Number of subjects
pars.K = 6; % Number of classes
pars.C = 3; % Number of channels

pars.runpar = 8; % The number of workers to use in parfor (if zero, uses just a regular for-loop)

pars.debuglevel = 2;
pars.figix      = 1;

pars.samp = 1.5; % Sampling size
pars.bs   = [1 1 1 1 1 1];

% Which parts of the algorithm to run
pars.do.w   	= 1;
pars.do.bf  	= 1;
pars.do.a0      = 0; 
pars.do.v0      = 0;
pars.do.pr      = 1;
pars.do.mu      = 1; pars.do.amu = 0;
pars.do.writemu = 0;

% Iteration numbers and stopping tolerance
pars.nitmain   = 50;
pars.nitcb     = 1;
pars.itstrtreg = [2 10 20]; % affine, small def., large def.
pars.tol       = 1e-4;
pars.rparam    = [0 0.005 1 0.25 1]*5; % Diffeomorphic regularisation

% pars.pthmu = path2spmtpm; % Use default SPM TPMs   
pars.pthmu = '';            % Generate TPMs from multiple subject brains

pars.alam = 1e-1;

pars.bflam  = 1e-3;
pars.bffwhm = 60;

% Preprocessing options----------------------------------------------------
pars.preproc.imload      = 0; % After the images have been preprocessed the first time, this option can be set to 1, to skip preprocessing the input images each time
pars.preproc.denoiseimg  = 0;
pars.preproc.cropimg     = 0;
pars.preproc.makenonneg  = 0;
pars.preproc.mnialign    = 0;
pars.preproc.resetorigin = 0;

run(pars);