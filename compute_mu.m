% TODO---------------------------------------------------------------------
% -Do not warp images beforehand
% -Push instead of interpolate
% -Extend to missing data (when using more than one modality) -> also in Kmeans
% -Regularise template update using spm_shoot_blur.m 
% -Add regularisation for W and n update in update_pr (for more stability)
% -Why is L decreasing if deg > 1
% -Start weights update later?

clear;

addpath('./core')
addpath('./util')

%--------------------------------------------------------------------------
% Initial parameters
%--------------------------------------------------------------------------

pars = [];

pars.imdir   = '/home/mbrud/Dropbox/PhD/Data/IXI-2D/';
pars.tempdir = './temp/'; % A temporary folder which will contain the preprocessed and warped images 

pars.N = 32; % Number of subjects
pars.K = 6; % Number of classes
pars.C = 3; % Number of channels

pars.runpar = 8; % The number of workers to use in parfor (if zero, uses just a regular for-loop)

pars.debuglevel = 3;
pars.figix      = 1;

pars.samp = 1.5; % Sampling size
pars.bs   = [1 1 1 1 1 1];

% Which parts of the algorithm to run
pars.do.w   	= 1;
pars.do.bf  	= 1;
pars.do.a0      = 0; % Doesn't work well?
pars.do.v0      = 1;
pars.do.pr      = 1;
pars.do.mu      = 1;
pars.do.writemu = 0;

% Iteration numbers and stopping tolerance
pars.nitmain    = 100; % 100
pars.nitcb      = 1;
pars.it2strtreg = [2 5 15]; % affine, small def., large def.
pars.tol        = 1e-4;
pars.rparam     = [0 0.005 1 0.25 1]; % Diffeomorphic regularisation

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