function [pth_template,uniform] = init_template(pth_template,V,K,dir_output,do_avg_template_dim,vx_tpm)
if nargin<6, vx_tpm = 1.5; end

if isempty(pth_template)
    uniform = true;
    if ~do_avg_template_dim
        % Use SPM's default template dimensions and orientation matrix
        %------------------------------------------------------------------
        
        pth_spm_template = fullfile(spm('dir'),'tpm','TPM.nii');  
        V_spm_template   = spm_vol(pth_spm_template);   
        d                = V_spm_template(1).dim;
        M                = V_spm_template(1).mat;
    else
        % Compute average dimensions and orientation matrix so that template will cover all images    
        %------------------------------------------------------------------
        
        % Collect orientation matrices and dimensions from all input images
        M    = numel(V);
        mats = [];
        dms  = [];    
        for m=1:M
            S = numel(V{m});
            N = numel(V{m}{1});

             for s=1:S
                for n=1:N
                    Nii = nifti(V{m}{s}(n).fname);
                    M   = Nii.mat;
                    mats = cat(3,mats,M);

                    d = size(Nii.dat);
                    dms = cat(1,dms,d);
                end
            end
        end

        % Compute average orientation matrix and dimensions
        [M,d] = compute_avg_mat(mats,dms);

        % Adjust template voxel size
        nvx = vx_tpm*ones(1,3);

        c = [1    1    1    1
             1    1    d(3) 1
             1    d(2) 1    1
             1    d(2) d(3) 1
             d(1) 1    1    1
             d(1) 1    d(3) 1
             d(1) d(2) 1    1
             d(1) d(2) d(3) 1]';

        tc = M(1:3,1:4)*c;
        if spm_flip_analyze_images, tc(1,:) = -tc(1,:); end

        mx = round(max(tc,[],2)');
        mn = round(min(tc,[],2)');

        M = spm_matrix(mn)*diag([nvx 1])*spm_matrix(-[1 1 1]);

        d = ceil((M\[mx 1]')');
        d = d(1:3);

        if spm_flip_analyze_images, M = diag([-1 1 1 1])*M; end
    end

    % Create and write TPMs to disk
    dir_template = fullfile(dir_output,'template');
    if exist(dir_template,'dir'), rmdir(dir_template,'s'); end; mkdir(dir_template);    

    pth_template  = fullfile(dir_template,'template.nii');   
    [pth,nam,ext] = fileparts(pth_template);

    img  = zeros(d,'single');
    vols = cell(K,1);
    for k=1:K    
        vols{k} = fullfile(pth,[nam num2str(k) ext]);
        create_nii(vols{k},img,M,16,'template');
    end
    clear img

    matlabbatch                       = cell(1,1);
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = pth_template;
    matlabbatch{1}.spm.util.cat.dtype = 16;
    spm_jobman('run',matlabbatch);

    delete(fullfile(pth,[nam '.mat']));

    for k=1:K
        delete(vols{k});
    end
else
    uniform = false;
end
%==========================================================================