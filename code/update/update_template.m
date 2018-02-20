function L = update_template(L,pth_template,obj,iter)

sched = exp(linspace(5,0,25));

Nii = nifti(pth_template);
vx  = vxsize(Nii.mat);
mu  = single(Nii.dat(:,:,:,:));

sparam = [0.01 0.02 1]; 
prm    = [vx, prod(vx)*[sparam(1) sched(min(iter,numel(sched)))*sparam(2:3)]];

[mu,L1] = spm_shoot_blur_wp_load(obj,mu,prm,iter);
L       = [L,L1];

if sum(~isfinite(mu(:))) || ~isfinite(L1), error('sum(~isfinite(logtpm(:))) || ~isfinite(dL)'); end

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = mu;     
%==========================================================================