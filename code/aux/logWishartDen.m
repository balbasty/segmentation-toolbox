function lZ = logWishartDen(W,nu)
% Log of normalising term of log W(W,nu)
M     = size(W,1);
if M >= nu+1,
   %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
    lZ = 0;
    return;
end
lGamM = M*(M-1)/4*log(pi);
for m=1:M,
    lGamM = lGamM+gammaln((nu+1-m)/2);
end
lZ    = 0.5*nu*(logdet(W) + M*log(2)) + lGamM;
%==========================================================================