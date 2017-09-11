function wf = warp(f,phi,bs)
if nargin<3, bs = [1 1 1 0 0 0]; end

d = size(phi);
d = d(1:3);

K = size(f,4);

wf = zeros([d K],'single');
for k=1:K
    wf(:,:,:,k) = spm_diffeo('bsplins',f(:,:,:,k),phi,bs);
end
wf(~isfinite(wf)) = 0;
%==========================================================================