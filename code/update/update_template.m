function L = update_template(L,pth_template,obj,iter)

s     = spm_shoot_defaults;
sched = s.sched;

Nii = nifti(pth_template);
vx  = vxsize(Nii.mat);
mu0 = single(Nii.dat(:,:,:,:));

sparam = [0 2 0]; 
scl    = sched(min(iter,numel(sched)));
prm    = [vx, prod(vx)*[sparam(1) scl*sparam(2) sparam(3)]];

[mu,L1] = spm_shoot_blur_wp_load(obj,mu0,prm,iter,true);
L       = [L,L1];

if 0
    [mu,L1] = spm_shoot_blur_wp(obj,mu0,prm,iter,12,true);
end

if sum(~isfinite(mu(:))) || ~isfinite(L1), error('sum(~isfinite(logtpm(:))) || ~isfinite(dL)'); end

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = mu;     
%==========================================================================