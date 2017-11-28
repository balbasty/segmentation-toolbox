function y = identity(d)
% Generate an identity transform of size d(1) x d(2) x d(3)
y = zeros([d(1:3) 3],'single');
[y(:,:,:,1),y(:,:,:,2),y(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
%=======================================================================