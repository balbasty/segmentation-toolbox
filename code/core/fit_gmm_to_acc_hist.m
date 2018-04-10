function obj = fit_gmm_to_acc_hist(obj,pars)
% Fit a GMM to an accumulative histogram of image intensities
% FORMAT obj = kmeans_on_hist(obj,pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
K = pars.K;
x = -2000:2000;
h = zeros(1,numel(x));
for m=1:numel(obj)
    if pars.dat{m}.segment.kmeans_hist
        S = numel(obj{m});
        for s=1:S
            img = single(obj{m}{s}.image(1).private.dat(:,:,:));
            msk = msk_modality(img,obj{m}{s}.modality,obj{m}{s}.trunc_ct);
            img = img(msk(:));
            
            h1 = hist(img(:),x);
            h  = h + h1;
        end
    end
end
h = h(900:end)';
x = x(900:end)';

[mg,mn,vr] = spm_imbasics('fit_gmm2hist',h,x,K);

for m=1:numel(obj)
    if pars.dat{m}.segment.kmeans_hist
        S = numel(obj{m});
        for s=1:S
            for k=1:K
                obj{m}{s}.segment.gmm.pr.m(:,k) = mn(k);
                obj{m}{s}.segment.gmm.pr.b(k) = mean(mg);
                obj{m}{s}.segment.gmm.pr.n(k) = mean(mg);
                obj{m}{s}.segment.gmm.pr.W(:,:,k) = 1./(vr(k)*mean(mg));

                obj{m}{s}.segment.gmm.po.m(:,k) = mn(k);
                obj{m}{s}.segment.gmm.po.b(k) = mean(mg);
                obj{m}{s}.segment.gmm.po.n(k) = mean(mg);
                obj{m}{s}.segment.gmm.po.W(:,:,k) = 1./(vr(k)*mean(mg));
            end
        end
    end
end
%==========================================================================