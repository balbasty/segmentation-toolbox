clear;

addpath(genpath('../code'))
addpath('/cherhome/mbrud/dev/distributed-computing')
addpath('/cherhome/mbrud/dev/auxiliary-functions')

%--------------------------------------------------------------------------
%% Set parameters
%--------------------------------------------------------------------------

pars.name         = 'seg-subj';
pars.dir_output   = '/data/mbrud/data-seg';
pars.dat          = {};

m = 1;
pars.dat{m}.dir_data = '/data/mbrud/images/CT/CHROMIS-preproc-rn-ss/';
pars.dat{m}.S = 1;
pars.dat{m}.modality = 'CT';
pars.dat{m}.segment.pth_prior = '/home/mbrud/Desktop/TEMP/processed-templates/mod-prior-CHROMIS-preproc-rn-ss.mat';
pars.dat{m}.segment.print_ll = true;
pars.dat{m}.segment.print_seg = true;
pars.dat{m}.segment.verbose = true;
pars.dat{m}.segment.do_write_res = true;
pars.dat{m}.segment.mrf = 1;
pars.dat{m}.segment.write_bf = [false true];
pars.dat{m}.segment.samp = 2;

pars.pth_template = '/home/mbrud/Desktop/TEMP/processed-templates/svmod-template.nii';
% spm_check_registration(pars.pth_template)

pars = pars_default(pars);

%--------------------------------------------------------------------------
%% Init
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

f1 = spm_select('FPList',obj{1}{1}.dir_write,'^m.*\.nii$');
f2 = spm_select('FPList',obj{1}{1}.dir_write,'^c.*\.nii$');
spm_check_registration(char({f1,f2}))