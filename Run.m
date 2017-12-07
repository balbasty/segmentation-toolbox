clear;

addpath(genpath('./code'))

%--------------------------------------------------------------------------
S = [Inf 100]; % Number of subjects
K = 16; % Number of classes (if a template is used, then K will be set to the number of classes in that template)

%--------------------------------------------------------------------------
% Options for running algorithm on the FIL cluster (Holly)
obj.run_on_holly     = true;
obj.holly.jnam       = 'seg-CR';
obj.holly.jnam_dummy = 'dummy-CR';

%--------------------------------------------------------------------------
% Folder for all algorithm data (not used when running on Holly)
obj.dir_data = '/data-scratch/mbrud/data/segmentation-toolbox-parfor';  

%--------------------------------------------------------------------------
% Define data cell array, which should contain the following:
% {'pth_to_images',num_subjects,'CT'/'MRI','healthy'/'non-healthy','pth_to_labels'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
im = {};

% 2D
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/CT-aged-2D-den',S(1),'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/CT-CHROMIS-2D-den',S(1),'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/CT-healthy-2D-den',S(2),'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/IXI-2D',S(2),'MRI','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/2D-Data/OASIS-long-2D',S,'MRI','healthy',''};

% 3D
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-aged-noneck-den',S,'CT','healthy',''};
im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck-den',S(1),'CT','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-healthy-noneck-den',S(2),'CT','healthy',''};
im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/IXI-noneck',S(2),'MRI','healthy',''};
% im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/OASIS-long-noneck',S,'MRI','healthy',''};

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
obj.preproc.do_preproc    = false; % Do preprocessing on input images
obj.preproc.is_DICOM      = false; % If input images are DICOM, converts DICOM to Nifti % TODO (work in progress)
obj.preproc.rem_corrupted = true; % Try to remove CT images that are corrupted (e.g. bone windowed)
obj.preproc.realign       = true; % Realign to MNI space
obj.preproc.crop          = true; % Remove data outside of head
obj.preproc.crop_neck     = true; % Remove neck (the spine, etc.)
obj.preproc.denoise       = true; % Denoise CT images

%--------------------------------------------------------------------------
% The distance (mm) between samples (for sub-sampling input data--improves speed)
obj.samp = 2;

%--------------------------------------------------------------------------
% Segmentation parameters
obj.lkp    = 1; % Number of gaussians per tissue
obj.vb     = false; % Use a variational Bayesian mixture model
obj.wp_reg = 1e0; % Bias weight updates towards 1

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
obj.nitermain = 200;
obj.tolmain   = 1e-4;
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
obj.vx_TPM   = obj.samp; % Voxel size of template to be estimated
obj.deg      = 2;        % Degree of interpolation when sampling template
obj.tiny     = 1e-4;     % Strength of Dirichlet prior used in template construction
obj.fwhm_TPM = 1e-2;     % Ad hoc smoothing of template (improves convergence)

%--------------------------------------------------------------------------
% For debugging
obj.verbose = 2;
obj.figix   = 1;

%--------------------------------------------------------------------------
% Folder for algorithm results (subfolder of dir_data)
obj.dir_res  = 'results';

%==========================================================================
%% Run algorithm
%==========================================================================

spm_preprocx_run(obj,im,K);