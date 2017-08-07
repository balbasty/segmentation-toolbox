function [Pf,mgm] = weighted_softmax(lnf,lnw,msk)
dm = size(lnf);
K  = dm(4);

if nargin<3, msk = ones([prod(dm(1:3)) 1],'logical'); end

Nmsk = nnz(msk);

lnf  = reshape(lnf,[prod(dm(1:3)) K]);
nlnf = zeros([Nmsk K]);
Pf   = zeros([Nmsk K]);
for k=1:K
    nlnf(:,k) = lnf(msk,k);
end

mgm = 0;
for k=1:K
    Pf(:,k) = exp(nlnf(:,k) + lnw(k));
    mgm     = mgm + Pf(:,k);
end

for k=1:K
    Pf(:,k) = Pf(:,k)./mgm; 
end

if sum(~isfinite(Pf(:)))
    warning('~isfinite(Pf(:))');
end

if nargout>1
    mgm = bsxfun(@rdivide,exp(nlnf),mgm);
    mgm = nansum(mgm,1);
end
%==========================================================================