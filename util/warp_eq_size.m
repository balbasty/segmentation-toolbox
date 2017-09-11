function [pthwf,mat,d] = warp_eq_size(pthf,samp,tempdir)

[N,C] = size(pthf);

% Compute average orientation matrix---------------------------------------
mat = [];
d   = [];
for n=1:N
    for c=1:C
        Nii = nifti(pthf{n,c});
        mat = cat(3,mat,Nii.mat);
        dm  = size(Nii.dat);
        if numel(dm)==2
           dm(3) = 1;
        end
        d = cat(1,d,dm);
    end
end

[mat,d] = compute_avg_mat(mat,d);

vx       = sqrt(sum(mat(1:3,1:3).^2));
st       = samp./vx;
F        = diag(st);
F(1:4,4) = 1;
mat      = mat*F;
d        = max(floor(d./st),1);

% Warp images--------------------------------------------------------------

folder = fullfile(tempdir,'warped');
if (exist(folder,'dir') == 0)
    mkdir(folder);
end

pthwf  = cell(N,C);
for n=1:N
    folder = fullfile(tempdir,'warped',['n' num2str(n)]);
    if (exist(folder,'dir') == 0)
        mkdir(folder);
    end
        
    for c=1:C
        Nii  = nifti(pthf{n,c});
        matf = Nii.mat;        
        f    = single(Nii.dat(:,:,:));             
      
        % Warp image-------------------------------------------------------
        phi = affine_transf(matf\mat,identity(d));
        f   = warp(f,phi);           

        % Mask-------------------------------------------------------------
        msk = get_msk(f);
        
        f(~msk) = 0;
        
        % Make images have simillar means----------------------------------                 
        a = 512/mean(f(:));
        f = f*a;
        
        % Write warped image-----------------------------------------------
        [~,nam,ext] = fileparts(Nii.dat.fname);

        pthwf{n,c} = fullfile(folder,['w' nam ext]);    
        create_nii(pthwf{n,c},f,mat,'float32','wf');
    end
end
%==========================================================================