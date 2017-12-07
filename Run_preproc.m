clear;

addpath(genpath('./code'))

%--------------------------------------------------------------------------
S = [Inf]; % Number of subjects

%--------------------------------------------------------------------------
% Folder for all algorithm data (not used when running on Holly)
obj.dir_data = '/data-scratch/mbrud/data/segmentation-toolbox-preproc';  

%--------------------------------------------------------------------------
% Define data cell array, which should contain the following:
% {'pth_to_images',num_subjects,'CT'/'MRI','healthy'/'non-healthy','pth_to_labels'}
%
% pth_to_images: Assumes that images are organised into per subject subfolders, i.e. if
% subject one has a T1 and a T2 image, then those images should be in a subfolder, for example, S1.
im = {};

im{end + 1} = {'/data-scratch/mbrud/images/Preprocessed/CT-CHROMIS-noneck',S(1),'CT','healthy',''};

browse_subjects(im{1}{1});

%--------------------------------------------------------------------------
% Run the algorithm in parallel by setting number of workers (Inf uses maximum number available)
obj.num_workers  = Inf;
obj.run_on_holly = false;

%--------------------------------------------------------------------------
% Preprocessing options
obj.preproc.do_preproc    = true; % Do preprocessing on input images
obj.preproc.is_DICOM      = false; % If input images are DICOM, converts DICOM to Nifti % TODO (work in progress)
obj.preproc.rem_corrupted = false; % Try to remove CT images that are corrupted (e.g. bone windowed)
obj.preproc.realign       = true; % Realign to MNI space
obj.preproc.crop          = true; % Remove data outside of head
obj.preproc.crop_neck     = true; % Remove neck (the spine, etc.)
obj.preproc.denoise       = false; % Denoise CT images

%%
V = load_and_process_images(obj,im);
manage_parpool(0);

%%
show_V(V)