function Lam = get_cpLam(cp)
[C,K] = size(cp.po.m);
Lam   = zeros(C,C,K);
for k=1:K
    Lam(:,:,k) = cp.po.n(k)*cp.po.W(:,:,k); % Mean of Wishart prior
end   
%==========================================================================