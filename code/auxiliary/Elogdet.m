function e = Elogdet(W,nu)
% E[log(det(Lambda))], where Lambda ~ W(W,nu)
N = size(W,1);
e = LogDet(W) + N*log(2);
for n=1:N
    e = e + psi((nu-n+1)/2);
end
%==========================================================================