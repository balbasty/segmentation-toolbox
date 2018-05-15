clear; close all;
 
VERBOSE = true;
 
%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
 
pth_distributed_toolbox      = '/home/mbrud/dev/WTCN-computational-anatomy-group/distributed-computing';
pth_auxiliary_functions      = '/home/mbrud/dev/WTCN-computational-anatomy-group/auxiliary-functions';
pth_batch_processing_toolbox = '/home/mbrud/dev/WTCN-computational-anatomy-group/batch-preprocessing-toolbox';
 
% addpath
addpath(genpath('./../code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)
addpath(genpath(pth_batch_processing_toolbox))
 
%--------------------------------------------------------------------------
%% First do some preprocessing
%--------------------------------------------------------------------------
 
pars = [];
 
pars.dat{1}.dir_data = '/data-scratch/mbrud/images/CT/Parashkev-CT-aged';
pars.dat{1}.dir_preproc = '/home/mbrud/Data/TESTING/segmentation-toolbox/img-delirium';
pars.dat{1}.modality = 'CT';
pars.dat{1}.S = 1;
% pars.dat{1}.preproc.do_realign2mni = true;
pars.dat{1}.preproc.do_crop = true;
pars.dat{1}.preproc.do_rem_neck = true;
% pars.dat{1}.preproc.do_denoise = true;
 
pars = preproc_default(pars);
 
pars = read_images_preproc(pars); 
obj0 = init_obj_preproc(pars);
 
fprintf('Preprocessing')
for s=1:numel(obj0{1})
    fprintf('.')
    process_subject_preproc(obj0{1}{s});
end
fprintf('\n')

%--------------------------------------------------------------------------
%% Then run segmentation
%--------------------------------------------------------------------------
 
pars = [];
 
pars.pth_template = '/home/mbrud/Data/TPM/CT-and-IXI/template.nii';
pars.dir_output = '/home/mbrud/Data/TESTING/segmentation-toolbox/seg-delirium';
pars.verbose = 0;
 
% pars.dat{1}.dir_data = obj0{1}{1}.dir_preproc;
% pars.dat{1}.modality = obj0{1}{1}.modality;

pars.dat{1}.dir_data = '/data/mbrud/images/CT/healthy_preproc-ra-cr-rn-vx/';
pars.dat{1}.modality = 'CT';
pars.dat{1}.S = 1;

pars.dat{1}.healthy = true;
pars.dat{1}.print_subj_info = VERBOSE;
% pars.dat{1}.S = 8;
 
pars.dat{1}.segment.pth_prior = '/home/mbrud/Data/TPM/CT-and-IXI/pr-CT.mat';
pars.dat{1}.segment.samp = 2;
pars.dat{1}.segment.mix_wp_reg = -1;
pars.dat{1}.segment.do_bf = false;
pars.dat{1}.segment.print_ll = VERBOSE;
pars.dat{1}.segment.verbose = VERBOSE;
pars.dat{1}.segment.niter = 12;
pars.dat{1}.segment.lkp.part = [1 1 1 1 2 2 3 4 5 6 7 8 9 9];
 
pars.dat{1}.write_res.do_write_res = true;
pars.dat{1}.write_res.mrf = 2;
pars.dat{1}.write_res.nmrf_its = 10;
 
pars = segment_default(pars); 
 
pars        = read_images_segment(pars); 
pars        = init_template(pars); 
[obj1,pars] = init_obj_segment(pars);
 
fprintf('Segmenting')
for s=1:numel(obj1{1})
    fprintf('.')
    process_subject_segment(obj1{1}{s},pars.fig);
end
fprintf('\n')

%--------------------------------------------------------------------------
%% Visualise results
%--------------------------------------------------------------------------
 
s = 3;
f1 = obj1{1}{s}.image(1).fname;
f2 = spm_select('FPList',obj1{1}{s}.dir_write,'^c.*\.nii$');
spm_check_registration(char({f1,f2}))