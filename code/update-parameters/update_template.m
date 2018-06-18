function L = update_template(L,obj,pars,iter)
% Update template (NIfTI stored in pth_template)
% FORMAT L = update_template(L,obj,pars,iter)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

pth_template = obj{1}{1}.pth_template;
sparam0      = pars.sparam0;
sparam1      = pars.sparam1;

Nii = nifti(pth_template);
vx  = spm_misc('vxsize',Nii.mat);
a   = single(Nii.dat(:,:,:,:));

sched = fliplr(sparam1:10:sparam0); 
reg   = sched(min(iter,numel(sched)));
prm   = [vx 0 reg reg];

if 1
    [a,L1] = spm_shoot_blur_wp_load(obj,a,prm,iter,1);
else
    [a,L1] = spm_shoot_blur_wp(obj,a,prm,iter,12,true);
end
L  = [L,L1];

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = a;     
%==========================================================================