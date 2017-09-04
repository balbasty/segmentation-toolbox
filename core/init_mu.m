function [lnmu,Mmu] = init_mu(pthf,pthmu,Mf,d,K,samp)

Mmu = Mf;

if isempty(pthmu)
    % Create uniform atlas
    lnmu = log(ones([d K],'single')/K);
else
    Kmu = numel(spm_vol(pthmu));
    
    if K~=Kmu
       error('K~=Kmu');
    end
           
    % Affinely register atlas from file to first image
    Affine = eye(4);
    V      = spm_vol(pthf{1,1});
    fwhm   = 0;    
    tpm    = spm_load_priors8(pthmu);
    
%     Affine = spm_maff8(V{1,1},samp,(fwhm+1)*16,tpm,Affine,'mni'); % Close to rigid
    
    M   = tpm.M\Affine*Mf;    
    phi = affine_transf(M,identity(d));
    phi = double(phi);
    x1  = phi(:,:,:,1);
    y1  = phi(:,:,:,2);
    z1  = phi(:,:,:,3);
    
    mutmp = spm_sample_priors8(tpm,x1,y1,z1);
    
    mu  = zeros([d K],'single');
    for k=1:K
       mu(:,:,:,k) = mutmp{k};
    end
end
%=======================================================================