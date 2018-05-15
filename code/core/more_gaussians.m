function gmm = more_gaussians(gmm,part)

% Adjust priors
gmm.pr = adjust(gmm.pr,part);

% Adjust posteriors
gmm.po = adjust(gmm.po,part);
%==========================================================================

%==========================================================================
function pars = adjust(pars,part)
N  = size(pars.m,1);
Kb = numel(pars.n);

m0 = reshape(pars.m,N,Kb);
n0 = pars.n;
b0 = pars.b;
W0 = reshape(pars.W,N,N,Kb);
 
m = zeros(size(m0));
n = zeros(size(n0));
b = zeros(size(b0));
W = zeros(size(W0));
  
for k=1:Kb
    kk = sum(part==k);
    w  = 1./(1+exp(-(kk - 1)*0.25)) - 0.5;
    
    vr0 = inv(n0(k)*W0(:,:,k));
    pr0 = inv(vr0*(1 - w));  
    
    b(part==k)     = b0(k);
    m(:,part==k)   = sqrtm(vr0)*randn(N,kk)*w + repmat(m0(:,k),[1,kk]);
    n(part==k)     = n0(k);
    W(:,:,part==k) = repmat(pr0/n0(k),[1 1 kk]);  
end
 
pars.m = m;
pars.n = n;
pars.b = b;
pars.W = W;
%==========================================================================