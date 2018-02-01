function [gmm,mg] = init_gmm(gmm,K,N,buf,obj,uniform,modality,ix_subj)
% Begin with moments:
if uniform, mom = spm_kmeans2mom(buf,K,modality,ix_subj);
else,       mom = compute_moments(buf,K);
end

if gmm.ml
    vr1     = zeros(N,N);
    [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm,vr1);  
else    
    gmm.pr  = obj.pr;     
    [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);  
end

mg = ones(K,1);
%======================================================================= 