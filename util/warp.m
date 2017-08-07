function y = warp(f,phi,ord)
if nargin<3, ord = [1 1 1 0 0 0]; end

f   = single(f);
phi = single(phi);

d   = size(phi);
y   = zeros(d(1:3),'single');
for k=1:size(f,4)
    f(:,:,:,k) = spm_diffeo('bsplinc',f(:,:,:,k),ord);
    y(:,:,:,k) = spm_diffeo('bsplins',f(:,:,:,k),phi,ord);
end
y(~isfinite(y)) = 0;
%==========================================================================