function L = update_template(L,pth_template,pth_obj,pth_H,pth_gr,iter)
M     = numel(pth_obj);
tot_S = 0;
for m=1:M
    S = numel(pth_obj{m});    
    for s=1:S, tot_S = tot_S + 1; end
end

shdef = spm_shoot_defaults;
smits = shdef.smits;
sched = shdef.sched;

Nii = nifti(pth_template);
vx  = vxsize(Nii.mat);
mu  = single(Nii.dat(:,:,:,:));
K   = size(mu,4);

% sparam_shoot = shdef.sparam;
% sparam       = sparam_shoot;
% prm          = [vx, prod(vx)*[sparam(1:2) sched(min(iter,numel(sched)))*sparam(3)]];
% prm(4)       = prm(4) + tot_S*K*1e-6;
sparam = [0.01 0.02 1]; 
prm    = [vx, prod(vx)*[sparam(1:2) sched(min(iter,numel(sched)))*sparam(3)]];
prm(4) = prm(4) + tot_S*K*1e-6;

% Update TPMs
if 0    
    t      = zeros([size(mu0) tot_S],'single');
    log_wp = zeros([tot_S K]);
    for m=1:M
        S = numel(pth_obj{m}); 
        for s=1:S
            obj          = load(pth_obj{m}{s}); 
            for k=1:K
                Nii          = nifti(obj.pth_resp{k});
                t(:,:,:,k,s) = single(Nii.dat(:,:,:));
            end
            log_wp(s,:)  = log(obj.wp);
        end
    end
    clear Nii

    [mu,~,dL] = spm_shoot_blur_wp(t,log_wp,prm,smits,mu0);
else
    [mu,~,dL] = spm_shoot_blur_wp_load(pth_obj,pth_H,pth_gr,mu,prm,iter);
end
L = [L,dL];

if sum(~isfinite(mu(:))) || ~isfinite(dL), error('sum(~isfinite(logtpm(:))) || ~isfinite(dL)'); end

% Save updated template
Nii              = nifti(pth_template);
Nii.dat(:,:,:,:) = mu;     
%==========================================================================