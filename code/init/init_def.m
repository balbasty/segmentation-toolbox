function [buf,param,MT,sk4,Twarp,llr] = init_def(buf,obj,lkp,sk,vx,ff,d,fig,wp,x0,y0,z0,tpm,M)
nz      = numel(buf);
Kb      = max(lkp);
pth_vel = obj.pth_vel;

spm_diffeo('boundary',1);

% This part is fiddly because of the regularisation of the warps.
% The fact that displacement fields are only parameterised every few
% voxels means that the functions in spm_diffeo need tweaking to
% account for the difference between the units of displacement and
% the separation of the voxels (if that makes sense).

% More work/thought is needed in terms of adjusting regularisation to
% account for different voxel sizes.  I'm still not satisfied that
% this (rescaling the regularisaiton by prod(vx.*sk)) is optimal.
% The same thing applies to all the nonlinear warping code in SPM.
param  = [sk.*vx prod(vx.*sk)*ff*obj.reg]; % FIX THIS (remove "prod(vx.*sk)")

% Mapping from indices of subsampled voxels to indices of voxels in image(s).
MT = [sk(1) 0 0 (1-sk(1));0 sk(2) 0 (1-sk(2)); 0 0 sk(3) (1-sk(3));0 0 0 1];

% For multiplying and dividing displacements to map from the subsampled voxel indices
% and the actual image voxel indices.
sk4 = reshape(sk,[1 1 1 3]);

if ~(exist(pth_vel,'file')==2)
    create_nii(pth_vel,zeros([d 3],'single'),eye(4),'float32','vel'); 
end

Nii   = nifti(pth_vel);
Twarp = zeros([d 3],'single');
for i=1:3
    Twarp(:,:,:,i) = single(Nii.dat(:,:,:,i)); 
end

llr = -0.5*sum(sum(sum(sum(Twarp.*bsxfun(@times,spm_diffeo('vel2mom',bsxfun(@times,Twarp,1./sk4),param),1./sk4)))));

% Load the warped prior probability images into the buffer
for z=1:nz
    if ~buf(z).Nm, continue; end
    
    msk1                          = buf(z).code>0;
    [x1,y1,z1]                    = make_deformation(Twarp,z,x0,y0,z0,M,msk1);
    b                             = spm_sample_priors8(tpm,x1,y1,z1);
    for k1=1:Kb, buf(z).dat(:,k1) = b{k1}; end
end

debug_view('template',fig{3},lkp,buf,wp);  
%=======================================================================