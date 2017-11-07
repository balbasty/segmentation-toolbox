function v = calc_var_of_grad(pth)
Nii        = nifti(pth);
img        = Nii.dat(:,:,:);
msk        = get_msk(img);
img(~msk)  = NaN;
vx         = sqrt(sum(Nii.mat(1:3,1:3).^2));
[Dx,Dy,Dz] = grad(img,vx);
gm         = sqrt(Dx(:).^2 + Dy(:).^2 + Dz(:).^2);

v = var(gm(isfinite(gm)));

% vf         = sum(img(isfinite(img)).^2) / sum(isfinite(img(:)));
         
% v         = (sum(Dx(isfinite(Dx)).^2) + sum(Dy(isfinite(Dy)).^2) + sum(Dz(isfinite(Dz)).^2)) ...
%             / (sum(isfinite(Dx(:))) + sum(isfinite(Dy(:))) + sum(isfinite(Dz(:))));         
               