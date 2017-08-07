function P = preprocess_im(imfname,tempdir,makenonneg,denoiseimg,cropimg,mnialign,resetorigin)  
if nargin<4, denoiseimg  = false; end
if nargin<5, cropimg     = false; end
if nargin<6, mnialign    = false; end
if nargin<7, resetorigin = false; end

N = numel(imfname);
P = cell(N,1);

f = fullfile(tempdir,'preproc');
if (exist(f,'dir') == 0)
    mkdir(f);
end

fprintf('Preprocessing image: '); 
for n=1:N
    fprintf('%d ',n); 
     
    P{n} = fullfile(f,['im' num2str(n) '.nii']);
   
    cpy_img(imfname{n},P{n});                                  
    
    if makenonneg
        nonneg(P{n},1000);
    end
        
    if resetorigin
        nm_reorient(P{n},1,1);
        reset_origin(P{n});
    end
    
    if mnialign
        rigid_align(P{n});    
    end
    
    if cropimg
        atlas_crop(P{n}); 
    end
    
    if denoiseimg
        denoise_img(P{n});            
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
