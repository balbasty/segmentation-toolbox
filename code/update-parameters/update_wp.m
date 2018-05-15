function wp = update_wp(lkp,s0,mgm,nvox,wp_reg,iter_template)
% Update tissue mixing weights (wp)
% FORMAT wp = update_wp(lkp,s0,mgm,nvox,wp_reg,iter_template)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
 
Kb = max(lkp.part);
wp = zeros(1,Kb);
for k1=1:Kb
    if strcmp(wp_reg,'single-subject')
        % Single subject segmentation -> default spm_preproc8 update
        wp(k1) = (sum(s0(lkp.part==k1)) + 1)/(mgm(k1) + Kb); % bias the solution towards 1
    elseif strcmp(wp_reg,'build-template')
        % Building template, decrease regularisation over iterations
        sched = 0:0.1:0.9;
        
        if iter_template>numel(sched)
            wp(k1) = (sum(s0(lkp.part==k1)) + 1)/(mgm(k1) + Kb);
        else
            w1 = sched(iter_template);
            w2 = 1 - w1;
 
            wp(k1) = (w1*sum(s0(lkp.part==k1)) + (w2*nvox)/Kb)/(w1*mgm(k1) + w2*nvox);
        end
    end
end
 
wp(lkp.rem) = eps;
wp          = wp/sum(wp);
%==========================================================================