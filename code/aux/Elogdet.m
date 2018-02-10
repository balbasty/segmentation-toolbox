function e = Elogdet(W,nu)
% E[log(det(Lambda))], where Lambda ~ W(W,nu)
M = size(W,1);
e = logdet(W) + M*log(2);
for m=1:M,
    e = e+psi((nu-m+1)/2);
end
%==========================================================================