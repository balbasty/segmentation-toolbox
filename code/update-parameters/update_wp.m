function wp = update_wp(lkp,s0,mgm,nvox,mix_wp_reg,wp_reg)
% Update tissue mixing weights (wp)
% FORMAT wp = update_wp(lkp,s0,mgm,nvox,mix_wp_reg,wp_reg)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Kb = max(lkp.part);
wp = zeros(1,Kb);
for k1=1:Kb
    if mix_wp_reg==-1
        % Default spm_preproc8 update
        wp(k1) = (sum(s0(lkp.part==k1)) + wp_reg*1)/(mgm(k1) + wp_reg*Kb); % bias the solution towards 1
    else
        % Weighted modfied version
        w1 = mix_wp_reg;
        w2 = 1 - w1;

        wp(k1) = (w1*sum(s0(lkp.part==k1)) + w2*nvox*wp_reg)/(w1*mgm(k1) + w2*nvox);
    end
end
wp = wp/sum(wp);
%==========================================================================