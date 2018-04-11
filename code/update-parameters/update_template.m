function L = update_template(L,obj,pars,iter)
pth_template = obj{1}{1}.pth_template;
sparam       = pars.sparam;

Nii = nifti(pth_template);
vx  = spm_misc('vxsize',Nii.mat);
mu0 = single(Nii.dat(:,:,:,:));

prm = [vx sparam];

if 1
    [mu,L1] = spm_shoot_blur_wp_load(obj,mu0,prm,iter,1);
else
    [mu,L1] = spm_shoot_blur_wp(obj,mu0,prm,iter,12,true);
end
L  = [L,L1];

if sum(~isfinite(mu(:))) || ~isfinite(L1), 
    error('sum(~isfinite(mu(:))) || ~isfinite(dL)'); 
end

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = mu;     
%==========================================================================