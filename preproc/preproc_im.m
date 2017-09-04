function pth = preproc_im(pthf,preproc,tempdir)  

[N,C] = size(pthf);
pth   = cell(N,C);

folder = fullfile(tempdir,'preproc');
if (exist(folder,'dir') == 0)
    mkdir(folder);
end

fprintf('Preprocessing image:\n'); 
for n=1:N
    fprintf('%d ',n);
    
    folder = fullfile(tempdir,'preproc',['n' num2str(n)]);
    if (exist(folder,'dir') == 0)
        mkdir(folder);
    end     
    
    for c=1:C        
        pth{n,c} = fullfile(folder,['im' num2str(c) '.nii']);

        cpy_img(pthf{n,c},pth{n,c});                                  

        if preproc.makenonneg
            nonneg(pth{n,c},1000);
        end

        if preproc.resetorigin
            nm_reorient(pth{n,c},1,1);
            reset_origin(pth{n,c});
        end

        if preproc.mnialign
            rigid_align(pth{n,c});    
        end

        if preproc.cropimg
            atlas_crop(pth{n,c}); 
        end

        if preproc.denoiseimg
            denoise_img(pth{n,c});            
        end
    end        
end
fprintf('\n'); 
%==========================================================================

%==========================================================================
function msk = nonneg(P,val)    

Nii = nifti(P);
img = Nii.dat(:,:,:);
dim = size(Nii.dat(:,:,:));    
mat = Nii.mat;

mn       = min(img(:));    
msk      = img == mn;
img(msk) = -val;        

img        = img + val;
V.fname    = P;
V.dim(1:3) = dim;
V.mat      = mat;
V          = spm_create_vol(V);
spm_write_vol(V,img);
%==========================================================================
