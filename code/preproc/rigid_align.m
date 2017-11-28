function R = rigid_align(P)
% Reposition an image by affine aligning to MNI space and Procrustes adjustment
% FORMAT rigid_align(P)
%        P - name of NIfTI image
% Image will have the matrix in its header adjusted.
% __________________________________________________________________________


% Load tissue probability data
tpm    = fullfile(spm('dir'),'tpm','TPM.nii,');
tpm    = [repmat(tpm,[6 1]) num2str((1:6)')];
tpm    = spm_load_priors8(tpm);

% Do the affine registration
Affine = eye(4);
Affine = spm_maff8(P,2,32,tpm,Affine,'mni'); % Heavily regularised
Affine = spm_maff8(P,2,1 ,tpm,Affine,'mni'); % Lightly regularised

% Load header
Nii    = nifti(P);

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(size(tpm.dat{1})),tpm.M);

% Generate mm coordinates of where deformation maps to
y1     = affind(x,inv(Affine));

% Weight the transform via GM+WM
weight = single(exp(tpm.dat{1})+exp(tpm.dat{2}));

% Weighted Procrustes analysis
[~,R]  = spm_get_closest_affine(x,y1,weight);

% Invert
% R      = inv(R);

% Write the new matrix to the header
Nii.mat = R\Nii.mat;
create(Nii);
%==========================================================================

%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
return;
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
return;
%==========================================================================