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
pars = '/home/mbrud/Dropbox/PhD/Data/pars/segmentation-toolbox/segment-CROMIS-3d.json';

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

f1 = spm_select('FPList',obj{1}{1}.dir_write,'^m.*\.nii$');
f2 = spm_select('FPList',obj{1}{1}.dir_write,'^c.*\.nii$');
spm_check_registration(char({f1,f2}))