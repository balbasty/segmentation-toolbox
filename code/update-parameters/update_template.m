function L = update_template(L,obj,pars,iter)
% Update template (NIfTI stored in pth_template)
% FORMAT L = update_template(L,obj,pars,iter)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

pth_template = obj{1}{1}.pth_template;
sparam       = pars.sparam;

Nii = nifti(pth_template);
vx  = spm_misc('vxsize',Nii.mat);
a   = single(Nii.dat(:,:,:,:));

sched       = fliplr(10:10:sparam); 
sparam(2:3) = sched(min(iter,numel(sched))); 
prm         = [vx sparam];

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