function [obj,K] = init_tpm(obj,V,K)
if isfield(obj,'pth_logtpm')
    pth_logtpm = obj.pth_logtpm;
else
    pth_logtpm     = '';
    obj.pth_logtpm = pth_logtpm;
end

if ~isempty(pth_logtpm) || numel(V{1})==1    
    % If not building TPMs (i.e., regular segmentation)
    obj.dopr        = false;
    obj.nitermain   = 1;    
    obj.num_workers = 0;    
    obj.dotpm       = false;    
    obj.niter       = 30;
    obj.niter1      = 8;
    obj.nsubitmog   = 20;
    obj.nsubitbf    = 1;
    obj.nitdef      = 3;
    obj.wp_reg      = 1;
end

if isempty(pth_logtpm)
    % Create uniform (log) TPMs
    [obj.pth_logtpm] = init_uniform_tpm(obj,V,K);
        
    obj.use_tpm = false;   
else    
    % TPMs will be loaded from file
    Nii     = nifti(pth_logtpm);
    d       = size(Nii.dat(:,:,:,:));
    clear Nii
    
    if (V{1}{1}(1).dim(3)>1 && d(3)==1) || (V{1}{1}(1).dim(3)==1 && d(3)>1)
        error('Template has incorrect dimensions!')
    end
    
    K           = d(4);    
    obj.use_tpm = true;         
end
%==========================================================================

%==========================================================================
function pth_logtpm = init_uniform_tpm(obj,V,K)
vxtpm      = obj.samp;
dir_res    = obj.dir_res;
pth_logtpm = fullfile(dir_res,'logtpm.nii');
       
% Create uniform TPMs that covers all images    
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

[pth,nam,ext] = fileparts(pth_logtpm);

img  = zeros(dm);
vols = cell(K,1);
for k=1:K    
    vols{k} = fullfile(pth,[nam num2str(k) ext]);
    create_nii(vols{k},img,mat,16,'logtpm');
end
clear img

matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_logtpm;
matlabbatch{1}.spm.util.cat.dtype = 16;
spm_jobman('run',matlabbatch);

delete(fullfile(pth,[nam '.mat']));

for k=1:K
    delete(vols{k});
end
%==========================================================================