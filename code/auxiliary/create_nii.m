function create_nii(fname,dat,mat,dtype,descrip)
Nii             = nifti;
d               = size(dat);
Nii.dat         = file_array(fname,d,dtype);
Nii.mat         = mat;
Nii.mat0        = mat;
Nii.mat_intent  = 'Scanner';
Nii.mat0_intent = 'Scanner';
Nii.descrip     = descrip;

create(Nii);

if numel(d)==4
    Nii.dat(:,:,:,:) = dat;
else
    Nii.dat(:,:,:) = dat;
end

