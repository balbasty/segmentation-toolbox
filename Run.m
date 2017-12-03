close all; clear;

addpath(genpath('./code'))

%--------------------------------------------------------------------------
S = 1; % Number of subjects
K = 6; % Number of classes (if a template is used, then K will be set to the number of classes in that template)

%--------------------------------------------------------------------------
% Options for running algorithm on the FIL cluster (Holly)
obj.run_on_holly     = 0;
obj.holly.jnam       = 'seg-h';
obj.holly.jnam_dummy = 'dummy-h';
obj.holly.RAM        = 6;
obj.holly.split      = 4;

%--------------------------------------------------------------------------
% Define data cell array, which should contain the following:
% {'pth_to_images',num_subjects,'CT'/'MRI','healthy'/'non-healthy','pth_to_labels'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
im = {};

% laptop
% im{end + 1} = {'/home/smajjk/Dropbox/PhD/Data/2D-Data/IXI-2D',S,'MRI','healthy',''};
im{end + 1} = {'/home/smajjk/Dropbox/PhD/Data/2D-Data/OASIS-Longitudinal-2D',S,'MRI','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/IXI-2D',S,'MRI','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/CT-CHROMIS-2D/',S,'CT','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/Rob-CT-healthy/',Inf,'CT','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/CT-aged-2D',S,'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/OASIS-Longitudinal-2D',S,'MRI','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-aged-noneck',S,'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/OASIS-noneck',S,'MRI','healthy',''};

% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck',S,'CT','healthy',''};

% browse_subjects(im{1}{1});

%--------------------------------------------------------------------------
% Path to initial template
obj.pth_logTPM = ''; % Path to existing template (set to '' for estimating a template, or get_spm_TPM for using the default SPM one)
% obj.pth_logTPM = '/home/smajjk/Dropbox/PhD/Data/logTPM/logTPM_2D.nii'; % Path to existing template (set to '' for estimating a template, or get_spm_TPM for using the default SPM one)

%--------------------------------------------------------------------------
% Run the algorithm in parallel by setting number of workers (Inf uses maximum number available)
obj.num_workers = Inf;

%--------------------------------------------------------------------------
% Preprocessing options
obj.preproc.do_preproc    = 0; % Do preprocessing on input images
obj.preproc.is_DICOM      = 0; % If input images are DICOM, converts DICOM to Nifti % TODO (work in progress)
obj.preproc.rem_corrupted = 1; % Try to remove CT images that are corrupted (e.g. bone windowed)
obj.preproc.realign       = 1; % Realign to MNI space
obj.preproc.crop          = 1; % Remove data outside of head
obj.preproc.crop_neck     = 1; % Remove neck (the spine, etc.)
obj.preproc.denoise       = 1; % Denoise CT images

%--------------------------------------------------------------------------
% The distance (mm) between samples (for sub-sampling input data--improves speed)
obj.samp = 2;

%--------------------------------------------------------------------------
% Segmentation parameters
obj.lkp    = 1; % Number of gaussians per tissue
obj.vb     = 1; % Use a variational Bayesian mixture model
obj.wp_reg = 1e0; % Bias weight updates towards 1

%--------------------------------------------------------------------------
% What estimates to perform
obj.dobias = 0; % Bias field
obj.doaff  = 1; % Affine registration
obj.dodef  = 1; % Non-linear registration
obj.dopr   = 1; % Intensity priors
obj.dotpm  = 1; % Template
obj.dowp   = 0; % Tissue mixing weights

%--------------------------------------------------------------------------
% Iterations and stopping tolerances for algorithm
obj.nitermain = 200;
obj.tolmain   = 1e-5;
obj.tolseg    = 1e-5;

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
obj.vx_TPM   = 1.5;  % Voxel size of template to be estimated
obj.deg      = 2;    % Degree of interpolation when sampling template
obj.tiny     = 1e-4; % Strength of Dirichlet prior used in template construction
obj.fwhm_TPM = 1e-2; % Ad hoc smoothing of template (improves convergence)
obj.mrf      = 0;    % Use a MRF cleanup procedure

%--------------------------------------------------------------------------
% For debugging
obj.verbose = 2;
obj.figix   = 1;

%--------------------------------------------------------------------------
% Make some folders
obj.dir_data = './data';  % for all algorithm data
% obj.dir_data = '/data-scratch/mbrud/data/segmentation-toolbox-parfor';  % for all algorithm data
obj.dir_res  = 'results'; % for algorithm results (subfolder of dir_data)

%==========================================================================
%% Run algorithm
%==========================================================================

spm_preprocx_run(obj,im,K);