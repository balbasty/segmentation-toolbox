function create_nii(pth,dat,mat,dtype,descrip)

if exist(pth,'file')==2, delete(pth); end

Nii         = nifti;
dm          = size(dat);
Nii.dat     = file_array(pth,dm,dtype);
Nii.mat     = mat;
Nii.mat0    = mat;
Nii.descrip = descrip;
create(Nii);

if numel(dm)==4
    Nii.dat(:,:,:,:) = dat;
else
    Nii.dat(:,:,:)   = dat;
end