function [Q,ll] = latent(f,bf,mg,gmm,B,lkp,wp,msk,code,labels,wp_lab,cr)
if nargin<12, cr = []; end

if isempty(labels)
    wp1 = 1;
    wp2 = 0;
else
    wp1 = 1 - wp_lab;
    wp2 = wp_lab;    
    
    tiny   = 1e-4;
    labels = log_spatial_priors(double(labels) + tiny,[],wp2);
end

B  = log_spatial_priors(B,wp,wp1);
Q  = log_likelihoods(f,bf,mg,gmm,msk,code,lkp,cr);

Kb   = max(lkp.part);
msk1 = code>0;
for k1=1:Kb
    for k=find(lkp.part==k1)
        if isempty(labels)
            Q(msk1,k) = Q(msk1,k) + B(:,k1);
        else
            Q(msk1,k) = Q(msk1,k) + B(:,k1) + labels(:,k1);
        end
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
        L(3) = nansum(nansum(Q(msk1,:).*bsxfun(@plus,B(:,lkp.part),log(mg)')));   
        
        ll   = sum(L);
    end
end
%=======================================================================