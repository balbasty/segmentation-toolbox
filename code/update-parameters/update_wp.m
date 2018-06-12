function wp = update_wp(lkp,s0,mgm,nvox,wp_reg)
% Update tissue mixing weights (wp)
% FORMAT wp = update_wp(lkp,s0,mgm,nvox,wp_reg)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
 
Kb = max(lkp.part);
wp = zeros(1,Kb) + eps;
for k1=1:Kb
    if strcmp(wp_reg,'single-subject')
        % Single subject segmentation -> default spm_preproc8 update
        wp(k1) = (sum(s0(lkp.part==k1)) + 1)/(mgm(k1) + Kb); % bias the solution towards 1
    elseif strcmp(wp_reg,'build-template')
        % Building template
        w1     = 0.8;
        w2     = 1 - w1;
        wp(k1) = (w1*sum(s0(lkp.part==k1)) + (w2*nvox)/Kb)/(w1*mgm(k1) + w2*nvox);
    end
end
 
wp(lkp.rem) = eps;
wp          = wp/sum(wp);
%==========================================================================