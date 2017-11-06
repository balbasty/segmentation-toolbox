function [Plogtpm,Kb,uniform,use_tpm] = init_template(Ptpm,V,Kb,vxtpm,dirTPM)
Plogtpm = fullfile(dirTPM,'TPM.nii');
if ~isempty(Ptpm)
    % Load template from file
    use_tpm = true;

    Vtpm    = spm_vol(Ptpm);    
    Kb      = numel(Vtpm);
    uniform = false;
    
    Nii = nifti(Ptpm);
    if V{1}{1}(1).dim(3)==1  
        % Subject data is 2D                
        d1  = floor(Vtpm(1).dim(3)/2 + 1);
        img = log(Nii.dat(:,:,d1,:));  % OBS: remove log once Ptpm contains logs
        mat = Vtpm(1).mat*spm_matrix([0 0 d1 - 1]);                              
    else
        % Subject data is 3D        
        img = log(Nii.dat(:,:,:,:)); % OBS: remove log once Ptpm contains logs
        mat = Nii.mat;
    end
    img(~isfinite(img)) = 0;
    clear Nii  
    
    vols = cell(Kb,1);
    for k=1:Kb
        [pth,nam,ext] = fileparts(Plogtpm);
        vols{k}       = fullfile(pth,[nam num2str(k) ext]);                        
        create_nii(vols{k},img(:,:,:,k),mat,Vtpm(1).dt(1),'TPM'); 
    end    
    clear img
    
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = Plogtpm;
    matlabbatch{1}.spm.util.cat.dtype = Vtpm(1).dt(1);
    spm_jobman('run',matlabbatch);
else       
    % No template provided so create a uniform TPM that covers all images    
    use_tpm = false;    
    uniform = true;

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
    vols = cell(Kb,1);
    for k=1:Kb
        [pth,nam,ext] = fileparts(Plogtpm);
        vols{k}       = fullfile(pth,[nam num2str(k) ext]);
        create_nii(vols{k},img,mat,16,'TPM');
    end
    clear img
    
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = Plogtpm;
    matlabbatch{1}.spm.util.cat.dtype = 16;
    spm_jobman('run',matlabbatch);
             
    [pth,nam] = fileparts(Plogtpm);
    delete(fullfile(pth,[nam '.mat']));
end

for k=1:Kb
    delete(vols{k});
end   
