function [po,pr] = vmstep(mom,pr)
if nargin<2, pr = []; end

mom = mom_John2Bishop(mom);  
s0  = mom.s0;
s1  = mom.s1;
S2  = mom.S2;

K = numel(s0);
N = size(s1,1);

if isempty(pr)
    m0 = s1;
    b0 = 0.01*ones(1,K);
    n0 = N*ones(1,K)-0.99;

    W0 = zeros(N,N,K);
    for k1=1:K
        W0(:,:,k1) = inv(S2(:,:,k1))/n0(k1)/1e+1;
    end
    
    pr.m = m0;
    pr.b = b0;
    pr.W = W0;
    pr.n = n0;
else
    m0 = pr.m;
    W0 = pr.W;
    b0 = pr.b;
    n0 = pr.n;
end

b = b0 + s0;
n = n0 + s0;
m = zeros(N,K);
W = zeros(N,N,K);
for k=1:K
    m(:,k) = (b0(k)*m0(:,k) + s0(k).*s1(:,k))./b(k);
    
    W0inv    = inv(W0(:,:,k));
    mlt1     = b0(k).*s0(k)/(b0(k) + s0(k));
    diff1    = s1(:,k) - m0(:,k);
    W(:,:,k) = inv(W0inv + s0(k)*S2(:,:,k) + mlt1*(diff1*diff1'));
    
    % Protect against eigenvalues becoming too "small"
    [V,D1]   = eig(W(:,:,k));
    tol      = max(diag(D1))*eps('single');
    D1       = diag(max(diag(D1),tol));
    W(:,:,k) = real(V*D1*V');
end  

po.m = m;
po.b = b;
po.W = W;
po.n = n;
%=======================================================================