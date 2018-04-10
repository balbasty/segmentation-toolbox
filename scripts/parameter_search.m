clear;

addpath(genpath('../../segmentation-toolbox'))

dir_output = '/data/mbrud/parameter-search/tmp-build-tpm-';
if exist(dir_output,'dir'), rmdir(dir_output,'s'); end; mkdir(dir_output);   

pars.verbose  = 0;                                                                                                                                                                   0;
pars.niter    = 30;
pars.dat{1}.S = 16;

vals = 1;

for i=1:numel(vals)
    spm_misc('manage_parpool',Inf);
          
    pars.dir_output = [dir_output num2str(i)];
    
%     pars.dat{1}.segment.mix_wp_reg = vals(i);
    
    build_template(pars,2);
    
    spm_misc('manage_parpool',0);
end

%% Visualise
i = 1;

f1 = figure(100);
f2 = figure(101);

pth_template = fullfile([dir_output num2str(i)],'res-test-local/template.nii');

show_template(f1,pth_template);

%%
pth_obj = fullfile([dir_output num2str(i)],'res-test-local/obj.mat');
obj     = load(pth_obj);
obj     = obj.obj;
show_resp(f2,obj,{1:8});