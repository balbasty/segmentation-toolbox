function [lnPx,ll] = latent(f,bf,mg,gmm,lnPa,lnPzN,lkp,wp,msk,code,labels,wp_l,cr)
if nargin<13, cr = []; end

Kb   = max(lkp.part); 
tiny = 1e-4; 
msk1 = code>0; 

lnPa = log_spatial_priors(lnPa,wp);
lnPx = log_likelihoods(f,bf,mg,gmm,msk,code,cr);

if ~isempty(labels) 
    lnPa1 = zeros([size(lnPa,1) Kb]); 
end

for k1=1:Kb
    for k=find(lkp.part==k1)
        if ~isempty(labels)
            if lkp.lab(k1)~=0               
                msk_l = labels==lkp.lab(k1);
 
                beta1 = log(1 - tiny); % log-probability of true labels
                beta2 = log(tiny);    % log-probability of false labels
 
                lnPa1(msk_l,k1)  = (1 - wp_l)*lnPa(msk_l,k1) + wp_l*beta1;
                lnPa1(~msk_l,k1) = (1 - wp_l)*lnPa(~msk_l,k1) + wp_l*beta2;                
            else
                lnPa1(:,k1) = lnPa(:,k1);
            end
            
            lnPx(msk1,k) = lnPx(msk1,k) + lnPa1(:,k1) + lnPzN(msk1,k1);
        else
            if numel(lnPzN)==Kb
                lnPx(msk1,k) = lnPx(msk1,k) + lnPa(:,k1) + lnPzN(:,k1);
            else
                lnPx(msk1,k) = lnPx(msk1,k) + lnPa(:,k1) + lnPzN(msk1,k1);
            end
        end
    end
end

if ~isempty(labels) 
    lnPa = lnPa1;
    clear B1
end

logSumQ = spm_matcomp('logsumexp',lnPx,2);
logQ    = bsxfun(@minus,lnPx,logSumQ);
lnPx    = exp(logQ);

if nargout==2
    % Compute lower bound components
    L = zeros(1,3);

    % Responsibilities
    L(1) = -nansum(nansum(lnPx.*logQ)) + nansum(nansum(bsxfun(@times,lnPx,log(mg)))); 

    % Bias field      
    N                        = numel(f);
    M                        = numel(code);
    nbf                      = ones([M N]);
    for n=1:N, nbf(msk{n},n) = double(bf{n}); end
    bf                       = nbf; clear nbf

    L(2) = nansum(nansum(bsxfun(@times,lnPx,log(prod(bf,2)))));

    % TPMs  
    L(3) = nansum(nansum(lnPx(msk1,:).*lnPa(:,lkp.part)));    
    
    ll = sum(L);
end
%=======================================================================