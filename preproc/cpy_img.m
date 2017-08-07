function cpy_img(P,nP)
Nii = nifti(P);
create_nii(nP,Nii.dat(:,:,:),Nii.mat,Nii.dat.dtype,Nii.descrip);
%==========================================================================