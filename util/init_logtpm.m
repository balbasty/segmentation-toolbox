function [Plogtpm,Kb] = init_logtpm(use_tpm,V,Kb,vxtpm,Plogtpm)
if nargin<5, Plogtpm = './tpm/TPM.nii'; end

Pspm    = which('spm');
Pspm    = Pspm(1:end-5);
Pspmtpm = fullfile(Pspm,'tpm','TPM.nii');
Vspmtpm = spm_vol(Pspmtpm);    

if use_tpm            
    Kb = numel(Vspmtpm);
    
    Nii = nifti(Pspmtpm);
    if V{1}(1).dim(3)==1  
        % Subject data is 2D                
        d1  = floor(Vspmtpm(1).dim(3)/2 + 1);
        img = log(Nii.dat(:,:,d1,:));        
        mat = Vspmtpm(1).mat*spm_matrix([0 0 d1 - 1]);                              
    else
        % Subject data is 3D        
        img = log(Nii.dat(:,:,:,:));
        mat = Nii.mat;
    end
    img(~isfinite(img)) = 0;
    clear Nii  
    
    vols = cell(Kb,1);
    for k=1:Kb
        [~,nam,ext] = fileparts(Plogtpm);
        fname       = [nam num2str(k) ext];                        
        create_nii(fname,img(:,:,:,k),mat,Vspmtpm(1).dt(1),'TPM');            
        vols{k}     = fname;
    end    
    clear img
    
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = Plogtpm;
    matlabbatch{1}.spm.util.cat.dtype = Vspmtpm(1).dt(1);
    spm_jobman('run',matlabbatch);

    for k=1:Kb
        delete(vols{k});
    end
else       
    % Create a uniform tpm that covers all images                
    S = numel(V);
    N = numel(V{1});
    
    mat = [];
    dm  = [];
    for s=1:S
        for c=1:N
            Nii = nifti(V{s}(c).fname);
            M   = Nii.mat;
            mat = cat(3,mat,M);
            
            d = size(Nii.dat);
            if numel(d)==2
               d(3) = 1;
            end        
            dm = cat(1,dm,d);
        end
    end

    [mat,dm] = compute_avg_mat(mat,dm,vxtpm);
    
    img  = log(ones(dm)/Kb);
    vols = cell(Kb,1);
    for k=1:Kb
        [~,nam,ext] = fileparts(Plogtpm);
        fname       = [nam num2str(k) ext];
        create_nii(fname,img,mat,Vspmtpm(1).dt(1),'TPM');
        vols{k} = fname;
    end
    clear img
    
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = Plogtpm;
    matlabbatch{1}.spm.util.cat.dtype = Vspmtpm(1).dt(1);
    spm_jobman('run',matlabbatch);

    for k=1:Kb
        delete(vols{k});
    end            
    
    [pth,nam] = fileparts(Plogtpm);
    delete(fullfile(pth,[nam '.mat']));
end