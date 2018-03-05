clear;

addpath(genpath('../code'))
addpath('/cherhome/mbrud/dev/distributed-computing')
addpath('/cherhome/mbrud/dev/auxiliary-functions')

%--------------------------------------------------------------------------
% Create image array
%--------------------------------------------------------------------------
S           = 1;
im          = {};       
im{end + 1} = {'/home/mbrud/Data/CT-big-lesions/',...
               S,'CT',[],3,'mean',...
               '/home/mbrud/Data/template/CT-lesion-local/prior-CT-big-lesions.mat',''};    
% im{end + 1} = {'/home/mbrud/Data/CT-big-lesions/',...
%                S,'CT',[],3,'random','',''};   
           
%--------------------------------------------------------------------------
%% Set parameters
%--------------------------------------------------------------------------
pars              = [];
pars.K            = 8;
pars.dir_output   = '/home/mbrud/Data/temp-data/seg_subj';
pars.dir_template = pars.dir_output;
pars.do_segment   = true;
pars.do_preproc   = false;
pars.vx_tpm       = 1.5;
pars.pth_template = '/home/mbrud/Data/template/CT-lesion-local/template.nii';
% pars.pth_template = '';

pars.segment.do_ml           = false;
pars.segment.do_bf           = true;
pars.segment.do_def          = true;
pars.segment.do_wp           = true;
pars.segment.do_write_res    = true;
pars.segment.do_push_resp    = false;
pars.segment.do_missing_data = false;
pars.segment.do_old_segment  = false;
pars.segment.kmeans_dist     = 'cityblock';
pars.segment.print_ll        = true;
pars.segment.print_seg       = true;    
pars.segment.verbose         = true;   
pars.segment.trunc_ct        = false;   
pars.segment.mrf             = 1;
% pars.segment.lkp             = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];

%--------------------------------------------------------------------------
%% Init
%--------------------------------------------------------------------------

V = cell(1,numel(im));
for m=1:numel(im) 
    V{m} = read_images(im{m},pars); 
end

pars = init_template(V,pars); 

[obj,pars,fig,rand_subjs] = init_obj(V,im,pars);

%--------------------------------------------------------------------------
%% Segment
%--------------------------------------------------------------------------

M = numel(obj);
for m=1:M
    S = numel(obj{m});
    for s=1:S
        update_subject(obj{m}{s},fig);
    end
end

%--------------------------------------------------------------------------
%% Visualise
%--------------------------------------------------------------------------

[files,dirs] = spm_select('FPList',obj{1}{1}.dir_write,'^c.*\.nii$');
spm_check_registration(files)