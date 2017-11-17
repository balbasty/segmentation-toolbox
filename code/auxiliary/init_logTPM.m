function pth_logTPM = init_logTPM(V,K,vxtpm,dir_res)
pth_logTPM = fullfile(dir_res,'logTPM.nii');
       
% Create a uniform TPM that covers all images    
M   = numel(V);
mat = [];
dm  = [];    
for m=1:M
    S = numel(V{m});
    N = numel(V{m}{1});

    for s=1:S
        for n=1:N
            Nii = nifti(V{m}{s}(n).fname);
            M   = Nii.mat;
            mat = cat(3,mat,M);

            d = size(Nii.dat);
            if numel(d)==2
               d(3) = 1;
            end        
            dm = cat(1,dm,d);
        end
    end
end

[mat,dm] = compute_avg_mat(mat,dm,vxtpm);

img  = zeros(dm);
vols = cell(K,1);
for k=1:K
    [pth,nam,ext] = fileparts(pth_logTPM);
    vols{k}       = fullfile(pth,[nam num2str(k) ext]);
    create_nii(vols{k},img,mat,16,'logTPM');
end
clear img

matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_logTPM;
matlabbatch{1}.spm.util.cat.dtype = 16;
spm_jobman('run',matlabbatch);

[pth,nam] = fileparts(pth_logTPM);
delete(fullfile(pth,[nam '.mat']));

for k=1:K
    delete(vols{k});
end