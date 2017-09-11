% TODO PRIO----------------------------------------------------------------
% -Double check equations in update_cp
% -Comment and clean update_cp
% -Speed up L calculation for different parts

% TODO FUTURE--------------------------------------------------------------
% -Push instead of interpolate
% -Change parameterisation of registration
% -Add regularisation for W and n update in update_pr (for more stability)
% -Run on MRI+CT
% -Include healthy CT
% -Regularise weight update somehow

% TEST---------------------------------------------------------------------
% -Try reg w/o weight updates to look if L decreases
% -cp iters: [1 2 3 4 5 6 7 8] 

% Qs-----------------------------------------------------------------------
% -Is Pf correct?
% -If for example VM end up in the wrong class initially, it has a hard
% time changing. Maybe this is the reason for the behaviour of the weights?
% -What is a good first chapter structure? i.e. what should be in
% background for example.

clear;

addpath('./core')
addpath('./util')

%--------------------------------------------------------------------------
% Initial parameters
%--------------------------------------------------------------------------

pars = [];

pars.imdir   = '/home/mbrud/Dropbox/PhD/Data/IXI-2D/';
pars.tempdir = '/home/mbrud/temp/'; % A temporary folder which will contain the preprocessed and warped images 

pars.N = 8; % Number of subjects
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
pars.do.a0      = 1; 
pars.do.amu     = 1;
pars.do.v0      = 1;
pars.do.pr      = 1;
pars.do.mu      = 1; 
pars.do.writemu = 0;

% Iteration numbers and stopping tolerance
pars.nitmain   = 100;
pars.nitcpbf   = 1;
pars.nitcp0    = [1 2 3 4 5 6 7 8];
pars.itstrtreg = [4 10 20]; % affine, small def., large def.
pars.tol       = 1e-4;

% pars.pthmu = path2spmtpm; % Use default SPM TPMs   
pars.pthmu = '';            % Generate TPMs from multiple subject brains

pars.bflam  = 1e-3;
pars.bffwhm = 60;

pars.alam   = 1e1;
pars.abasis = 12;

pars.vlam = [0 0.005 1 0.25 1]; % Diffeomorphic regularisation

pars.mufwhm = 0.25;

% Preprocessing options----------------------------------------------------
pars.preproc.imload      = 1; % After the images have been preprocessed the first time, this option can be set to 1, to skip preprocessing the input images each time
pars.preproc.denoiseimg  = 0;
pars.preproc.cropimg     = 0;
pars.preproc.makenonneg  = 0;
pars.preproc.mnialign    = 0;
pars.preproc.resetorigin = 0;

run(pars);