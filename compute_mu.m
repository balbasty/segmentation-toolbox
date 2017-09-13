% TODO PRIO----------------------------------------------------------------
% -Double check equations in update_cp
% -Speed up L calculation for different parts
% -Use NaN instead of zeros for missing data

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
% -BF-correction for CT not improving estimates?
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

imdir{1}   = '/home/mbrud/Data/CT-w-lesion/';
imdir{2}   = '/home/mbrud/Dropbox/PhD/Data/CT-w-lesion-2D/';
imdir{3}   = '/home/mbrud/Dropbox/PhD/Data/IXI-2D/';
imdir{4}   = '/home/mbrud/Data/IXI-subjects/';
imdir{5}   = '/home/mbrud/Dropbox/PhD/Data/CT-healthy-2D/';
pars.imdir = imdir{5};

pars.tempdir = '/home/mbrud/temp/MRI/'; % A temporary folder which will contain the preprocessed and warped images 

pars.N = 100; % Number of subjects
pars.K = 10; % Number of classes
pars.C = 3;  % Number of channels

pars.runpar = 1; 

pars.ct = 1;

pars.preproc.imload = 0;

pars.debuglevel = 4;
pars.figix      = 1;

pars.samp = 1.5; % Sampling size
pars.bs   = [1 1 1 1 1 1];

% Which parts of the algorithm to run
pars.do.w   	= 1;
pars.do.bf  	= 0;
pars.do.a0      = 1; 
pars.do.amu     = 0;
pars.do.v0      = 1;
pars.do.pr      = 1;
pars.do.mu      = 1; 
pars.do.writemu = 0;

% Iteration numbers and stopping tolerance
pars.nitmain    = 200;
pars.nitcpbf    = 1;
pars.nitcp0     = [1 2 3 4 5 6 7 8];
pars.nitsrtw    = 1;
pars.nitstrtreg = [2 10 15]; % affine, small def., large def.
pars.tol        = 1e-4;

% pars.pthmu = path2spmtpm; % Use default SPM TPMs   
pars.pthmu = '';            % Generate TPMs from multiple subject brains

pars.bflam  = 1e-3;
pars.bffwhm = 60;

pars.alam   = 1e1;
pars.abasis = 12;

pars.vlam = [0 0.005 1 0.25 1]; % Diffeomorphic regularisation

pars.mufwhm = 0.2;

% Preprocessing options----------------------------------------------------
pars.preproc.mnialign    = 0; % Always do?
pars.preproc.denoiseimg  = 1;
pars.preproc.cropimg     = 0;
pars.preproc.resetorigin = 0; % Always do?

run(pars);