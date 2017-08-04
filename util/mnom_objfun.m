function ll = mnom_objfun(r,lnmu,y1,lnw,fig,ord,lnwmu)
if nargin<7, lnwmu  = []; end

d  = size(r);
ll = 0;

for z=1:d(3)
    ll = ll + mnom_objfun_slice(r,lnmu,y1,z,lnw,fig,ord,lnwmu);
end
%==========================================================================