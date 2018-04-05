function L = update_template(L,obj,pars,iter)
pth_template = obj{1}{1}.pth_template;
sparam       = pars.sparam;
s            = spm_shoot_defaults;
sched        = s.sched;
sched        = sched(2:end);
sched        = repelem(sched,2);

Nii = nifti(pth_template);
vx  = vxsize(Nii.mat);
mu0 = single(Nii.dat(:,:,:,:));

scl = sched(min(iter,numel(sched)));
% scl = 1;
prm = [vx, scl*prod(vx)*[sparam(1) sparam(2) sparam(3)]];

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