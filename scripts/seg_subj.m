clear;

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
% pth_distributed_toolbox = '/home/mbrud/dev/distributed-computing';
% pth_auxiliary_functions = '/home/mbrud/dev/auxiliary-functions';
pth_distributed_toolbox = '/cherhome/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/cherhome/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------

addpath(genpath('./../code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------
pars.pth_template = '/home/mbrud/Dropbox/PhD/Data/TPMs/template.nii';
pars.dir_output = '/data/mbrud/images/TESTING/OUTPUT/default';

pars.dat{1}.dir_data = '/data/mbrud/images/TESTING/CROMIS';
pars.dat{1}.modality = 'CT';
pars.dat{1}.print_subj_info = true;
pars.dat{1}.S = 1;

pars.dat{1}.segment.pth_prior = '/home/mbrud/Dropbox/PhD/Data/TPMs/prior-scans.mat';
pars.dat{1}.segment.samp = 3;
pars.dat{1}.segment.mix_wp_reg = -1;
pars.dat{1}.segment.do_bf = false;
pars.dat{1}.segment.print_ll = true;
pars.dat{1}.segment.verbose = true;
pars.dat{1}.segment.lkp.part = [1 1 1 1 1 1 1 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 10 10 10 10 10 10 10 10];

pars.dat{1}.write_res.do_write_res = true;
pars.dat{1}.write_res.mrf = 0;

pars = segment_default(pars); 

%--------------------------------------------------------------------------
% Initialise algorithm
%--------------------------------------------------------------------------
pars       = read_images(pars); 
pars       = init_template(pars); 
[obj,pars] = init_obj(pars);

%--------------------------------------------------------------------------
%% Segment
%--------------------------------------------------------------------------

process_subject(obj{1}{1},pars.fig);

%--------------------------------------------------------------------------
%% Visualise
%--------------------------------------------------------------------------

f1 = obj{1}{1}.image(1).fname;
f2 = spm_select('FPList',obj{1}{1}.dir_write,'^c.*\.nii$');
spm_check_registration(char({f1,f2}))