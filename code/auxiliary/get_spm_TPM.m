function P = get_spm_TPM
Pspm = which('spm');
Pspm = Pspm(1:end-5);
P    = fullfile(Pspm,'tpm','TPM.nii');