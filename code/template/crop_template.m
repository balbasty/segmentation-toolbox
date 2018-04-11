function crop_template(pth_template,iter,verbose)
% Crop template to same size as default SPM template
% FORMAT crop_template(pth_template,iter,verbose)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin<3, verbose = true; end

pth0 = fullfile(spm('dir'),'tpm','TPM.nii');  
V0   = spm_vol(pth0);   
d0   = V0(1).dim;
vx0  = spm_misc('vxsize',V0(1).mat);

V1  = spm_vol(pth_template);
K1  = numel(V1);
d1  = V1(1).dim;
vx1 = spm_misc('vxsize',V1(1).mat);

msk     = d1<d0;
d0(msk) = d1(msk);

sk0 = d0;
sk1 = d1;

bb1 = floor((sk1 - sk0)/2);
bb2 = bb1 + sk0;
bb  = [bb1' bb2'];

for k=1:K1
    spm_impreproc('subvol',V1(k),bb','');        
end

if verbose
    fprintf('%2d | size(otpm) = [%d %d %d] | size(ntpm) = [%d %d %d]\n',iter,d1(1),d1(2),d1(3),d0(1),d0(2),d0(3));
end
%==========================================================================