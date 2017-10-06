function V = reg_and_reslice(V)
S = numel(V);
D = numel(V{1});

matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';

if D>1
    for s=1:S
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {V{s}(1).fname};

        for d=2:D
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {V{s}(d).fname};
            
            delete(V{s}(d).fname);
            
            output_list = spm_jobman('run',matlabbatch);
            V{s}(d)     = spm_vol(output_list{1}.rfiles{1});
        end
    end
end
%==========================================================================