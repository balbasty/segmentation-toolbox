function [Q,ll] = latent(f,bf,mg,gmm,B,lkp,wp,msk,code,K_lab,cr)
if nargin<11, cr = []; end

B  = log_spatial_priors(B,wp);
Q  = log_likelihoods(f,bf,mg,gmm,msk,code,K_lab,lkp,cr);

Kb   = max(lkp);
msk1 = code>0;
for k1=1:Kb
    for k=find(lkp==k1)
        Q(msk1,k) = Q(msk1,k) + B(:,k1);
    end
end

if gmm.ml
    [Q,ll] = safe_softmax(Q);
else
    logSumQ = logsumexp(Q,2);
    logQ    = bsxfun(@minus,Q,logSumQ);
    Q       = exp(logQ);
    
    if nargout==2
        % Compute lower bound components
        L = zeros(1,3);

        % Responsibilities
        L(1) = -nansum(nansum(Q.*logQ));

        % Bias field      
        N                        = numel(f);
        M                        = numel(code);
        nbf                      = ones([M N]);
        for n=1:N, nbf(msk{n},n) = double(bf{n}); end
        bf                       = nbf; clear nbf
     
        L(2) = nansum(nansum(bsxfun(@times,Q,log(prod(bf,2)))));

        % TPMs
        L(3) = nansum(nansum(Q(msk1,:).*bsxfun(@plus,B(:,lkp),log(mg)')));   
        
        ll   = sum(L);
    end
end
%=======================================================================