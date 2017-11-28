function denoise_img(fname,descrip)

% Get input data
Nii = nifti(fname);
mat = Nii.mat;
vx  = sqrt(sum(mat(1:3,1:3).^2));
Y   = single(Nii.dat(:,:,:));   
dm  = size(Y);
if numel(dm)==2
   dm(3) = 1; 
end
clear Nii

% Define projection matrix
pm.A   = @(x) x;
pm.At  = @(x) x;
pm.AtA = @(x) x;

% Regularisation
lambda = 4; 

% Options
opts          = [];
opts.verbose  = 0;
opts.usegpu   = false;
opts.nonneg   = false;
opts.vx       = vx;
opts.bregmitr = 1;
opts.rho      = 0.25;
opts.lam2     = 1e-3;
opts.maxitr   = 100;
opts.reltol   = 1e-4;
opts.lstol    = 1e-4;

% Get mask values
msk  = get_msk(Y,descrip);
vals = single(Y(~msk));

Xhat = ADMM_img_solver(pm,{Y},lambda,opts);
clear Y

% Remask because values previously zero may have been slightly modified
Xhat{1}(~msk) = vals;
clear msk vals

% Write image
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['den_' nam ext]);
                
V.fname    = nfname;
V.dim(1:3) = dm;
V.mat      = mat;
V.dt       = [4 0];
V          = spm_create_vol(V);
spm_write_vol(V,Xhat{1});
%==========================================================================