function [buf,param,MT,sk4,Twarp,llr] = init_def_and_dat(buf,obj,sk,vx,ff,d,fig,wp,x0,y0,z0,tpm,M)
nz      = numel(buf);
Kb      = max(obj.segment.lkp.part);
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
param  = [sk.*vx prod(vx.*sk)*ff*obj.segment.reg]; % FIX THIS (remove "prod(vx.*sk)")

% Mapping from indices of subsampled voxels to indices of voxels in image(s).
MT = [sk(1) 0 0 (1-sk(1));0 sk(2) 0 (1-sk(2)); 0 0 sk(3) (1-sk(3));0 0 0 1];

% For multiplying and dividing displacements to map from the subsampled voxel indices
% and the actual image voxel indices.
sk4 = reshape(sk,[1 1 1 3]);

if ~(exist(pth_vel,'file')==2)
    spm_misc('create_nii',pth_vel,zeros([d 3],'single'),eye(4),obj.dt,'vel');
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
    b                             = spm_sample_logpriors(tpm,x1,y1,z1);
    for k1=1:Kb, buf(z).dat(:,k1) = b{k1}; end
end

debug_view('template',fig{3},obj.segment.lkp,buf,wp);  
%==========================================================================

%==========================================================================
function [x1,y1,z1] = make_deformation(Twarp,z,x0,y0,z0,M,msk)
x1a = x0    + double(Twarp(:,:,z,1));
y1a = y0    + double(Twarp(:,:,z,2));
z1a = z0(z) + double(Twarp(:,:,z,3));
if nargin>=7
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
if numel(z0)==1
   z1 = ones(size(z1));
end
return;
%==========================================================================