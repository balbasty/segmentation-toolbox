function mrf_template(pth_template,verbose)
if nargin<2, verbose = false; end

Nii = nifti(pth_template);
mat = Nii.mat;
vx  = 1./single(sum(mat(1:3,1:3).^2));
Q   = single(Nii.dat(:,:,:,:));
dm  = size(Q);
Kb  = dm(4);
zix = floor(dm(3)/2) + 1;

% softmax
Q = exp(Q);
Q = bsxfun(@rdivide,Q,sum(Q,4));

nmrf_its = 10;
T        = 2;
G        = T*ones([Kb,1],'single');
P        = zeros(dm,'uint8');

if verbose
    figure(666);
    for k=1:Kb
       subplot(2,Kb,k) 
       imagesc(Q(:,:,zix,k)); axis off image xy; colormap(gray);
    end
end

for iter=1:nmrf_its
    spm_mrf(P,Q,G,vx);
end

P = double(P)/255;

if verbose
    for k=1:Kb
       subplot(2,Kb,Kb + k) 
       imagesc(P(:,:,zix,k)); axis off image xy; colormap(gray);
    end
    drawnow
end

Nii.dat(:,:,:,:) = log(max(P,eps('single')));
%==========================================================================