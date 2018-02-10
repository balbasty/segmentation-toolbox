function L = update_template(L,pth_template,obj,iter)

sched = exp(linspace(5,0,25));

Nii = nifti(pth_template);
vx  = vxsize(Nii.mat);
mu  = single(Nii.dat(:,:,:,:));

sparam = [0.01 0.02 1]; 
prm    = [vx, prod(vx)*[sparam(1) sparam(2) sched(min(iter,numel(sched)))*sparam(3)]];

[mu,~,dL] = spm_shoot_blur_wp_load(obj,mu,prm,iter);
L         = [L,dL];
    
if sum(~isfinite(mu(:))) || ~isfinite(dL), error('sum(~isfinite(logtpm(:))) || ~isfinite(dL)'); end

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = mu;     
%==========================================================================