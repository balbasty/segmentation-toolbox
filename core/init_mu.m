function [mu,matmu] = init_mu(Px,Pmu,mat,d,K,samp)

if isempty(Pmu)
    % Create flat atlas
    mu    = ones([d K],'single')/K;
    matmu = mat;
else
    Kmu = numel(spm_vol(Pmu));
    
    if K~=Kmu
       error('K~=Kmu');
    end
           
    % Affinely register atlas from file to first image
    V      = spm_vol(Px{1});
    fwhm   = 0;
    Affine = eye(4);
    tpm    = spm_load_priors8(Pmu);
    affreg = 'mni';
    
    Affine = spm_maff8(V,samp,(fwhm+1)*16,tpm,Affine,affreg); % Close to rigid
    
    [x0,y0,z0] = ndgrid(1:d(1),1:d(2),1:d(3));    
    matmu      = mat;
    M          = tpm.M\Affine*mat;
    Twarp      = zeros([d,3],'single');
    [x1,y1,z1] = defs(Twarp,x0,y0,z0,M);
    mutmp      = spm_sample_priors8(tpm,x1,y1,z1);
    
    mu  = zeros([d K],'single');
    for k=1:K
       mu(:,:,:,k) = mutmp{k};
    end
end
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(Twarp,x0,y0,z0,M,msk)
x1a = x0 + double(Twarp(:,:,:,1));
y1a = y0 + double(Twarp(:,:,:,2));
z1a = z0 + double(Twarp(:,:,:,3));
if nargin>=7
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================