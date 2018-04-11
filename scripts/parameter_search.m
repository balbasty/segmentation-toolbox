clear;

pth2_distributed_toolbox = '/home/mbrud/dev/distributed-computing';
pth2_auxiliary_functions = '/home/mbrud/dev/auxiliary-functions';

addpath(pth2_distributed_toolbox)
addpath(pth2_auxiliary_functions)
addpath(genpath('../../segmentation-toolbox'))

dir_output = '/home/mbrud/Data/temp-segmentation-toolbox';
if exist(dir_output,'dir'), rmdir(dir_output,'s'); end; mkdir(dir_output);   

pars.verbose  = 0;                                                                                                                                                                   0;
pars.niter    = 30;
pars.dat{1}.S = 8;

vals = 0:0.1:1;
par  = 'mix_wp_reg'; % CT-0.9

pars.dat{1}.segment.do_bf = false;

for i=1:numel(vals)
    spm_misc('manage_parpool',Inf);
          
    pars.dir_output = [dir_output '-' par '-' num2str(i)];
    
    pars.dat{1}.segment.mix_wp_reg = vals(i);
    
    build_template(pars,2);
    
    spm_misc('manage_parpool',0);
end

vals = 10.^-(0:6);
par  = 'bf';

pars.dat{1}.segment.do_bf = true;
pars.dat{1}.segment.mix_wp_reg = 0.5;

for i=1:numel(vals)
    spm_misc('manage_parpool',Inf);
          
    pars.dir_output = [dir_output '-' par '-' num2str(i)];
    
    pars.dat{1}.segment.biasreg = vals(i);
    
    build_template(pars,2);
    
    spm_misc('manage_parpool',0);
end

%% Visualise
i   = 10;
par = 'mix_wp_reg';

pth_template = fullfile([dir_output '-' par '-' num2str(i)],'res-test-local/template.nii');

f1 = figure(100);
show_template(f1,pth_template);

%%
pth_obj = fullfile([dir_output '-' par '-' num2str(i)],'res-test-local/obj.mat');
obj     = load(pth_obj);
obj     = obj.obj;

f2 = figure(101);
show_resp(f2,obj,{1:8});