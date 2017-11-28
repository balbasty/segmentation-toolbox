function V = reg_and_reslice(V)
N = numel(V);

% Get image with largest volume and reslice using this image as reference
vol = zeros(N,3);
for n=1:N
    vx       = sqrt(sum(V(n).mat(1:3,1:3).^2));
    vol(n,:) = vx.*V(n).dim;
end
vol        = prod(vol,2);
[~,ref_ix] = max(vol);

% Set options
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 2;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {V(ref_ix).fname};

% Register and reslice
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {V(n).fname};                        

    output_list = spm_jobman('run',matlabbatch);
    delete(V(n).fname);

    V(n) = spm_vol(output_list{1}.rfiles{1});            
end
%==========================================================================