function [Q,ll] = latent(f,bf,mg,gmm,B,lkp,wp,code,cr)
if nargin<9, cr = []; end

B  = log_spatial_priors(B,wp);
Q  = log_likelihoods(f,bf,mg,gmm,code,cr);
Kb = max(lkp);
for k1=1:Kb
    for k=find(lkp==k1)
        Q(:,k) = Q(:,k) + B(:,k1);
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
        L(1) = -sum(sum(Q.*logQ));

        % Bias field    
        N                   = numel(f);
        M                   = numel(f{1});
        nbf                 = ones([M N]);
        for n=1:N, nbf(:,n) = double(bf{n}); end
        bf                  = nbf; clear nbf
     
        L(2) = sum(sum(bsxfun(@times,Q,log(prod(bf,2)))));

        % TPMs
        L(3) = sum(sum(Q.*bsxfun(@plus,B(:,lkp),log(mg)')));   
        
        ll   = sum(L);
    end
end
%=======================================================================