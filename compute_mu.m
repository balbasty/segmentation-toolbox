% TODO---------------------------------------------------------------------
% -Correct parameterisation of registration parameters?
% -Do not warp images to same size
% -Push instead of interpolate
% -Extend to missing data (when using more than one modality)
% -Make compatible with multi-channel data
% -Does affine registration really work well?
% -Correct precision for bias field estimation?

% Qs-----------------------------------------------------------------------
% -If multiple classes per tissue, how to chose number of classes?
% -For tissue updates, OK to skip update if nL<oL

addpath('./core')
addpath('./util')

%--------------------------------------------------------------------------
% Initial parameters
%--------------------------------------------------------------------------

pars = [];

pars.N = 8; % Number of subjects
pars.K = 6;  % Number of classes

pars.imdir   = '/home/smajjk/Data/IXI-T1/'; % A folder containing nifti MRIs
pars.tempdir = './tempmri/';                         % A temporary folder which will contain the preprocessed and warped images 

% Preprocessing options----------------------------------------------------
pars.imload      = 0; % After the images have been preprocessed the first time, this option can be set to 1, to skip preprocessing the input images each time
pars.denoiseimg  = 0;
pars.cropimg     = 0;
pars.makenonneg  = 0;
pars.mnialign    = 0;
pars.resetorigin = 0;

% pars.pthmu = path2spmtpm; % Use default SPM TPMs   
pars.pthmu = '';            % Generate TPMs from multiple subject brains

pars.samp = 3; % Sampling size
pars.ord  = 1; % Degree of b-spline interpolation

% Set to 1 to display TPMs and some random bias fields and velocity fields
pars.debugmode  = 1;
pars.plotL      = 1;
pars.showwarp   = 1;

% Which parts of the algorithm to run
pars.docp    = 1;
pars.dobias  = 1;
pars.doaff   = 1;
pars.doprior = 1;
pars.dotpm   = 1;

% Iteration numbers and stopping tolerance
pars.nitmain    = 50;
pars.nitcb      = 8;
pars.nitaff     = 4;
pars.startvelit = 2;
pars.tol        = 1e-4;

pars.int_args = 8;                  % Set to 1 for small deformation approximation
pars.rparam   = [0 0.005 1 0.25 1]; % Diffeomorphic regularisation

run(pars);