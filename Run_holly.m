clear;

addpath(genpath('./code'))

%--------------------------------------------------------------------------
S = [72 72]; % Number of subjects
K = 16; % Number of classes (if a template is used, then K will be set to the number of classes in that template)

%--------------------------------------------------------------------------
% Options for running algorithm on the FIL cluster (Holly)
obj.run_on_holly = true;
obj.holly_jnam   = 'CT-MRI';
obj.holly_ram    = 6;

%--------------------------------------------------------------------------
% Define data cell array, which should contain the following:
% {'pth_to_images',num_subjects,'CT'/'MRI','healthy'/'non-healthy','pth_to_labels'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
im = {};

im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-healthy-neck',S(1),'CT','healthy',''};
im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/IXI-neck',S(2),'MRI','healthy',''};

%--------------------------------------------------------------------------
% Run the algorithm in parallel by setting number of workers (Inf uses maximum number available)
obj.num_workers = Inf;

%--------------------------------------------------------------------------
% Segmentation parameters
obj.samp         = 2; % The distance (mm) between samples (for sub-sampling input data--improves speed)
obj.missing_data = true;
obj.vb           = true; % Use a variational Bayesian mixture model
obj.wp_reg       = 1e0; % Bias weight updates towards 1
obj.lkp          = 1; % Number of gaussians per tissue

%--------------------------------------------------------------------------
% What estimates to perform
obj.dobias = true; % Bias field
obj.doaff  = true; % Affine registration
obj.dodef  = true; % Non-linear registration
obj.dopr   = true; % Intensity priors
obj.dotpm  = true; % Template
obj.dowp   = true; % Tissue mixing weights

%--------------------------------------------------------------------------
% Iterations and stopping tolerances for algorithm
obj.nitermain = 100;
obj.tolmain   = 1e-5;
obj.tolseg    = 1e-4;

obj.niter     = 2;
obj.niter1    = 20;
obj.nsubitmog = 20;
obj.nsubitbf  = 1;
obj.nitdef    = 3;

%--------------------------------------------------------------------------
% Regularisation for estimating deformations
obj.rparam = [0 0.001 0.5 0.05 0.2]*0.1;

%--------------------------------------------------------------------------
% Bias field options
obj.biasfwhm = 60;
obj.biasreg  = 1e-4;       

%--------------------------------------------------------------------------
% Template options
obj.deg          = 2;    % Degree of interpolation when sampling template
obj.pr_dirichlet = 1e-4; % Strength of Dirichlet prior used in template construction
obj.fwhmtpm      = 1e-2; % Ad hoc smoothing of template (improves convergence)
obj.crop_bb      = true;

%--------------------------------------------------------------------------
% For debugging
obj.verbose = 2;
obj.distfig = false;
obj.figix   = 1;

%--------------------------------------------------------------------------
% Folder for algorithm results (subfolder of dir_data)
obj.dir_res  = 'results';

%--------------------------------------------------------------------------
% Preprocessing options
obj.preproc.do_preproc    = false; % Do preprocessing on input images
obj.preproc.rem_corrupted = false; % Try to remove CT images that are corrupted (e.g. bone windowed)

%==========================================================================
%% Run algorithm
%==========================================================================
spm_preprocx_run(obj,im,K);