function lZ = logWishartDen(W,nu)
% Log of normalising term of log W(W,nu)
N     = size(W,1);
if N >= nu+1
   %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
    lZ = 0;
    return;
end
lGamM = N*(N-1)/4*log(pi);
for n=1:N
    lGamM = lGamM + gammaln((nu+1-n)/2);
end
lZ = 0.5*nu*(LogDet(W) + N*log(2)) + lGamM;
%==========================================================================