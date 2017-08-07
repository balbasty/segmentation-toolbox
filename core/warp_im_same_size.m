function [pthx,pthmsk,mat,d] = warp_im_same_size(Px,samp,tempdir)

% Create folders, if doens not exist---------------------------------------
fw = fullfile(tempdir,'warped');
if (exist(fw,'dir') == 0)
    mkdir(fw);
end

fmsk = fullfile(tempdir,'msk');
if (exist(fmsk,'dir') == 0)
    mkdir(fmsk);
end

% Compute average orientation matrix---------------------------------------
N   = numel(Px);
mat = [];
d   = [];
for n=1:N
    Nii = nifti(Px{n});
    mat = cat(3,mat,Nii.mat);
    d   = cat(1,d,size(Nii.dat));
end

[mat,d] = compute_avg_mat(mat,d);

vx       = sqrt(sum(mat(1:3,1:3).^2));
st       = samp./vx;
F        = diag(st);
F(1:4,4) = 1;
mat      = mat*F;
d        = floor(d./st);

% Warp images and create masks---------------------------------------------
pthx   = cell(N,1);
pthmsk = cell(N,1);
for n=1:N 
    Nii  = nifti(Px{n});
    matn = Nii.mat;
        
    x = Nii.dat(:,:,:);             
    
    % Warp image-----------------------------------------------------------
    phi = affine_transf(matn\mat,identity(d));
    x   = warp(x,phi);
    
    % Create mask----------------------------------------------------------
    msk = x~=0 & x~=max(x(:)) & x~=min(x(:)) & isfinite(x) & x~=-3024 & x~=-1500;    
    
    % Make images have simillar means--------------------------------------                       
    a = 1024/mean(x(msk));
    x = x*a;
    
    % Write warped image and mask------------------------------------------
    [~,nam,ext] = fileparts(Nii.dat.fname);
    
    pthx{n} = fullfile(fw,['w' nam ext]);    
    create_nii(pthx{n},x,mat,'float32','X');
    
    pthmsk{n} = fullfile(fmsk,['msk' nam ext]);
    create_nii(pthmsk{n},msk,mat,'uint8-le','Mask');
end
%==========================================================================