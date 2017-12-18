function [Q,L] = latent(buf,B,mg,mog,wp,lkp,cr)
if nargin<7, cr = []; end
K = numel(lkp);
B = log_spatial_priors(B,wp,mg,lkp);
Q = log_likelihoods(K,buf,[],mog,cr);
Q = Q + B;

if isfield(mog,'pr')            
    logSumQ = logsumexp(Q,2);
    logQ    = bsxfun(@minus,Q,logSumQ);
    Q       = exp(logQ);
    
    % Compute lower bound components
    L = zeros(1,3);
    
    % Responsibilities
    L(1) = -nansum(nansum(Q.*logQ));
    
    % Bias field    
    N   = numel(buf.f);
    nbf = ones([numel(buf.msk{1}) N]);
    for n=1:N
        nbf(buf.msk{n},n) = double(buf.bf{n});
    end
    bf   = nbf; clear nbf    
    bf   = repmat(log(prod(bf,2)),1,K);
    L(2) = nansum(nansum(Q.*bf));
    
    % TPMs
    L(3) = sum(sum(Q.*B));   
else
    [Q,L] = safe_softmax(Q);
end
%=======================================================================

%=======================================================================
function B = log_spatial_priors(B,wp,mg,lkp)
%     if Kb~=K      
% %         mu  = bsxfun(@times,mu(:,lkp),mg);        
% %         for k=1:Kb
% %             kk  = find(lkp==k);
% %             smu = sum(mu(:,kk),2);
% %             for k1=kk                
% %                 mu(:,k1) = mu(:,k1)./smu;            
% %             end
% %         end
% %         mu    = bsxfun(@times,mu,wp(lkp)); 
% %         logmu = log(mu);
%         mu    = bsxfun(@times,mu,wp);
%         mu    = bsxfun(@times,mu,1./sum(mu,2));
%         mu    = bsxfun(@times,mu(:,lkp),mg);
%         logmu = log(mu);
%     end

B = bsxfun(@times,B,wp);
B = log(bsxfun(@times,B,1./sum(B,2)));
B = B(:,lkp);
B = bsxfun(@plus,B,log(mg));
%=======================================================================