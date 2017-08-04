function denoise_img(fname,lambda)
if nargin < 2, lambda = 2; end

% Get input data
Nii = nifti(fname);    
vx  = sqrt(sum(Nii.mat(1:3,1:3).^2));
Y   = Nii.dat(:,:,:);   
        
% Define projection matrix
pm.A   = @(x) x;
pm.At  = @(x) x;
pm.AtA = @(x) x;

% Set denoising options
opts          = [];
opts.verbose  = 0;
opts.usegpu   = false;
opts.nonneg   = false;
opts.vx       = vx;
opts.bregmitr = 1;
opts.rho      = 1e1;
opts.lam2     = 1e-3;
opts.maxitr   = 50;
opts.reltol   = 1e-4;
opts.lstol    = 1e-4;

% Get zero values
msk = Y == 0;

Xhat = ADMM_img_solver(pm,{Y},lambda,opts);

% Remask because values previously zero may have been slightly modified
Xhat{1}(msk) = 0;

% Write image
V.fname    = fname;
V.dim(1:3) = size(Y);
V.mat      = Nii.mat;
V          = spm_create_vol(V);
spm_write_vol(V,Xhat{1});
%==========================================================================