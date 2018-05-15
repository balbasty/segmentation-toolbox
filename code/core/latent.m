function [Q,ll] = latent(f,bf,mg,gmm,B,lkp,wp,msk,code,labels,wp_l,cr)
if nargin<12, cr = []; end

Kb   = max(lkp.part); 
tiny = 1e-4; 
msk1 = code>0; 

B = log_spatial_priors(B,wp);
Q = log_likelihoods(f,bf,mg,gmm,msk,code,cr);

if ~isempty(labels) 
    B1 = zeros([size(B,1) Kb]); 
end

for k1=1:Kb
    for k=find(lkp.part==k1)
        if ~isempty(labels)
            if lkp.lab(k1)~=0               
                msk_l = labels==lkp.lab(k1);
 
                beta1 = log(1 - tiny); % log-probability of true labels
                beta2 = log(tiny);    % log-probability of false labels
 
                B1(msk_l,k1)  = (1 - wp_l)*B(msk_l,k1) + wp_l*beta1;
                B1(~msk_l,k1) = (1 - wp_l)*B(~msk_l,k1) + wp_l*beta2;                
            else
                B1(:,k1) = B(:,k1);
            end
            
            Q(msk1,k) = Q(msk1,k) + B1(:,k1);
        else
            Q(msk1,k) = Q(msk1,k) + B(:,k1);
        end
    end
end

if ~isempty(labels) 
    B = B1;
    clear B1
end

logSumQ = spm_matcomp('logsumexp',Q,2);
logQ    = bsxfun(@minus,Q,logSumQ);
Q       = exp(logQ);

if nargout==2
    % Compute lower bound components
    L = zeros(1,3);

    % Responsibilities
    L(1) = -nansum(nansum(Q.*logQ)) + nansum(nansum(bsxfun(@times,Q,log(mg)))); 

    % Bias field      
    N                        = numel(f);
    M                        = numel(code);
    nbf                      = ones([M N]);
    for n=1:N, nbf(msk{n},n) = double(bf{n}); end
    bf                       = nbf; clear nbf

    L(2) = nansum(nansum(bsxfun(@times,Q,log(prod(bf,2)))));

    % TPMs  
    L(3) = nansum(nansum(Q(msk1,:).*B(:,lkp.part)));    

    ll = sum(L);
end
%=======================================================================