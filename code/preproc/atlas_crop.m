function Affine = atlas_crop(P,Affine,prefix,crop_neck)
if nargin < 2, Affine    = []; end
if nargin < 3, prefix    = ''; end
if nargin < 4, crop_neck = false; end

% Locate TPM.nii in SPM
pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');

Vin  = spm_vol(P);
Vtpm = spm_vol(pth_tpm);

mat    = Vin.mat;
mattpm = Vtpm.mat;

if isempty(Affine)
    tpm    = spm_load_priors8(Vtpm);
    samp   = 3;
    Affine = spm_maff8(P,samp,16,tpm,eye(4),'mni');
    Affine = spm_maff8(P,samp,16,tpm,Affine,'mni');
end

% Voxel locations in TPM.nii
Ltpm1 = [120 72.2 37.3 1]'; Ltpm2 = [120 72.2 75.9 1]';
Rtpm1 = [3  72.2 37.3 1]'; Rtpm2 = [3  72.2 75.9 1]';

Stpm1 = [58.6 42.6 119 1]'; Stpm2 = [58.60 99.4 119 1]';
if crop_neck
    Itpm1 = [58.6 39.4 2.5   1]'; Itpm2 = [58.60 99.4 2.5  1]';    
else
    Itpm1 = [58.6 39.4 -200  1]'; Itpm2 = [58.60 99.4 -200 1]';
end
Atpm1 = [58.6 144 28.4 1]'; Atpm2 = [58.60 144 82.3 1]';
Ptpm1 = [58.6 3.5  28.4 1]'; Ptpm2 = [58.60 3.5 82.3 1]'; 

% Voxel locations in input
T  = mat\(Affine\mattpm);
L1 = T*Ltpm1; L2 = T*Ltpm2;
R1 = T*Rtpm1; R2 = T*Rtpm2;
U1 = T*Stpm1; U2 = T*Stpm2;
D1 = T*Itpm1; D2 = T*Itpm2;
A1 = T*Atpm1; A2 = T*Atpm2;
P1 = T*Ptpm1; P2 = T*Ptpm2;

% Bounding-box
for i=1:3
    X     = [L1(i) R1(i) U1(i) D1(i) A1(i) P1(i)...
             L2(i) R2(i) U2(i) D2(i) A2(i) P2(i)];
    mx(i) = max(X);
    mn(i) = min(X);
end

% Do cropping
subvol(Vin,[mn(1) mx(1);mn(2) mx(2);mn(3) mx(3)]',prefix);
%==========================================================================