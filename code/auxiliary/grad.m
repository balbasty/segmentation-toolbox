function [Dx,Dy,Dz] = grad(X,vx) 
if nargin < 2, vx = ones(3,1); end

precision = get_type(X);

if size(X,3)==1
    Dx = [diff(X,1,2),zeros(size(X,1),1,precision)]./vx(2);
    Dy = [diff(X,1,1);zeros(1,size(X,2),precision)]./vx(1);
    Dz = 0;
else
    Dx = cat(2,diff(X,1,2),zeros(size(X,1),1,size(X,3),precision))./vx(2);
    Dy = cat(1,diff(X,1,1),zeros(1,size(X,2),size(X,3),precision))./vx(1);
    Dz = cat(3,diff(X,1,3),zeros(size(X,1),size(X,2),1,precision))./vx(3);  
end
%==========================================================================
  
%==========================================================================
function out = get_type(var)
tmp = whos('var');
out = tmp.class;
%==========================================================================