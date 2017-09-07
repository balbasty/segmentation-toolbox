function [phi,Jphi,theta,Jtheta] = make_deformation(v,prm,int_args,Greens)
% phi    - Forward deformation field n1*n2*n3*3 (single prec. float)
% Jphi   - Forward Jacobian tensors n1*n2*n3 (single prec. float)
% theta  - Inverse deformation field n1*n2*n3*3 (single prec. float)
% Jtheta - Inverse Jacobian tensors n1*n2*n3 (single prec. float)
if ~isempty(Greens)
    [phi,Jphi,~,theta,Jtheta] = spm_shoot3d(v,prm,int_args,Greens);
else
    dm = size(v);
    
    theta        = cell(3,1);
    [theta{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
    theta        = cat(4,theta{:});
    
    phi = theta + v;
    
    id                = ones(dm(1:3),'single');    
    Jtheta            = zeros([dm 3],'single');
    Jtheta(:,:,:,1,1) = id;
    Jtheta(:,:,:,2,2) = id;
    Jtheta(:,:,:,3,3) = id;
    
    Jphi = Jtheta;
end
%==========================================================================