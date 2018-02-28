function pars = init_template(V,K,pars)

if exist(pars.dir_template,'dir'), rmdir(pars.dir_template,'s'); end; mkdir(pars.dir_template);    

if ~pars.do_segment
    pars.pth_template = fullfile(spm('dir'),'tpm','TPM.nii');  
    V1                = spm_vol(pars.pth_template); 
    pars.uniform      = false;
    pars.dt           = V1.dt;
    return
end

if isempty(pars.pth_template)
    pars.uniform = true;
    
    % Compute average dimensions and orientation matrix so that template 
    % will cover all images    
    %----------------------------------------------------------------------

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
    nvx = pars.vx_tpm*ones(1,3);

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

    pars.dt = [spm_type('float32') spm_platform('bigend')];
    
    pars.pth_template = fullfile(pars.dir_template,'template.nii');   
    [pth,nam,ext]     = fileparts(pars.pth_template);

    img  = zeros(d,'single');
    vols = cell(K,1);
    for k=1:K    
        vols{k} = fullfile(pth,[nam num2str(k) ext]);
        create_nii(vols{k},img,M,pars.dt,'template');
    end
    clear img

    matlabbatch                       = cell(1,1);
    matlabbatch{1}.spm.util.cat.vols  = vols;
    matlabbatch{1}.spm.util.cat.name  = pars.pth_template;    
    matlabbatch{1}.spm.util.cat.dtype = 0;
    spm_jobman('run',matlabbatch);

    delete(fullfile(pth,[nam '.mat']));

    for k=1:K
        delete(vols{k});
    end
else
    pars.uniform = false;
end
%==========================================================================