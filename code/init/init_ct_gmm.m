function pars = init_ct_gmm(pars)
% Init GMM when data is CT
% FORMAT pars = init_ct_gmm(pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

Kb = pars.K; 
M  = numel(pars.dat);

% Read CT populations into new struct (pars_ct)
tot_S   = 0;
pars_ct = {};
cnt     = 1;
for m=1:M
    if isfield(pars.dat{m},'modality')
        modality = pars.dat{m}.modality;
    else
        continue
    end
    if strcmp(modality,'CT')
        pars_ct.dat{cnt} = pars.dat{m};
        tot_S            = pars_ct.dat{cnt}.S;
        cnt              = cnt + 1;
    end
end

if ~isempty(pars_ct) && tot_S>1
    % There are data that is CT -> perform CT init routine (below)
    %----------------------------------------------------------------------
    
    % Calculate histogram for each CT image and store the histograms
    c   = {};
    x   = {};
    cnt = 1;
    mn  = Inf;
    mx  = -Inf;
    for m=1:M
        S = pars_ct.dat{m}.S;
        for s=1:S
            Nii = nifti(pars_ct.dat{m}.V{s}.fname);
            img = single(Nii.dat(:,:,:));
            msk = msk_modality(img,'CT');
            img = img(msk(:));

            x1 = min(img(:)):max(img(:));
            c1 = hist(img(:),x1);

            if max(x1)>mx
                mx = max(x1);
            end

            if min(x1)<mn
                mn = min(x1);
            end

            c{cnt} = single(c1);
            x{cnt} = single(x1);

            cnt = cnt + 1;
        end
    end
    clear c1 x1

    % Sum all histograms into (H)
    X   = mn:mx;
    C   = zeros(1,numel(X));
    msk = false(1,numel(X));
    for i=1:numel(c)
        x1      = x{i};
        x1      = x1 + abs(mn) + 1;
        msk(x1) = true;

        C(msk) = C(msk) + c{i};

        msk(1,:) = false;
    end
    clear c x

    % Smooth histogram a little bit
    krn = spm_smoothkern(4,(-256:256)',0);
    C   = conv(C,krn,'same');
    
    % Remove intensity values smaller than a threshold (nmn)
    nmn = -1040;
    nmx = 2000;
    nix = find(X>=nmn & X<=nmx);
    X   = X(nix);
    C   = C(nix);

    % Fit VB-GMM to background
    [~,ix]          = find(X<0);
    [~,m1,b1,n1,W1] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),6);      
    lkp1            = ones(1,6);
    
    % Fit VB-GMM to brain
    k               = Kb - 2;
    [~,ix]          = find(X>=0 & X<=80);
    [~,m2,b2,n2,W2] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k); 
            
    ngauss = pars_ct.dat{m}.segment.ct.ngauss; % Number of Gaussians for each brain tissue class
    if ngauss>1
        % Use more than one Gaussian per class within the brain. This
        % requires correcting the Gauss-Wishart parameters estimated by the
        % histogram fit.
        vr = 1./(n2.*W2); % variance
                
        nm2 = [];
        nb2 = [];
        nn2 = [];
        nW2 = [];
        
        w  = 1./(1+exp(-(ngauss - 1)*0.25)) - 0.5;
        pr = 1./(vr*(1 - w));  
        
        for k1=1:k                        
            m22 = sqrtm(vr(k1))*randn(1,ngauss)*w + repmat(m2(k1),[1,ngauss]);   
            nm2 = [nm2,m22];
            
            b22 = b2(k1)/ngauss*ones(1,ngauss);
            nb2 = [nb2,b22];
            
            n22 = n2(k1)/ngauss*ones(1,ngauss);
            nn2 = [nn2,n22];
                        
            W22 = pr(k1)/n2(k1)*ones(1,ngauss);
            nW2 = [nW2,W22];
        end

        m2 = nm2;
        b2 = nb2;
        n2 = nn2;
        W2 = nW2;
    end
    lkp2 = repelem(2:Kb - 1,1,ngauss);
    
    % Fit VB-GMM to bone
    [~,ix]          = find(X>80);
    [~,m3,b3,n3,W3] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),6);     
    lkp3            = Kb*ones(1,6);
     
    m    = [m1,m2,m3];
    W    = [W1,W2,W3];
    b    = [b1,b2,b3];
    n    = [n1,n2,n3];
    part = [lkp1,lkp2,lkp3];
    K    = numel(part);
    
    [~,ix_rem] = min(abs(m - 65)); % Get index of lesion class
    rem        = part(ix_rem);
  
    % Init VB-GMM posteriors
    %----------------------------------------------------------------------
    W   = reshape(W,1,1,K);
    m   = reshape(m,1,K);
    gmm = struct;
    mg  = zeros(1,K);
    for m1=1:M
        S = numel(pars_ct.dat{m1});
        for s=1:S
            for k=1:Kb
                ix = find(part==k);

                mg(ix) = 1/nnz(ix);

                gmm.pr.m(:,ix)   = m(ix);
                gmm.pr.b(ix)     = 1;%b(ix);
                gmm.pr.n(ix)     = 1;%n(ix);
                gmm.pr.W(:,:,ix) = 1;%W(ix);

                gmm.po.m(:,ix)   = m(ix);
                gmm.po.b(ix)     = 1;%b(ix);
                gmm.po.n(ix)     = 1;%n(ix);
                gmm.po.W(:,:,ix) = 1;%W(ix);                
            end
        end
    end

    % Set CT GMM structs
    %----------------------------------------------------------------------    
    for m=1:M
        modality = pars.dat{m}.modality;
        healthy  = pars.dat{m}.healthy;
        if strcmp(modality,'CT')
            S = numel(pars_ct.dat{m});
            for s=1:S
                pars.dat{m}.segment.mg       = mg;
                pars.dat{m}.segment.lkp.part = part;                
                pars.dat{m}.segment.gmm      = gmm;
                
                if healthy
                    % If labelleled healthy, remove class with intensity
                    % closest to intensity of blood in CT. Also set that
                    % class' posterior and prior m hyper-parameter to zero.
                    pars.dat{m}.segment.lkp.rem  = rem;
                    
                    for i=1:ngauss
                        pars.dat{m}.segment.gmm.po.m(ix_rem + i - 1) = 0;
                        pars.dat{m}.segment.gmm.pr.m(ix_rem + i - 1) = 0;
                    end
                end
            end
        end
    end
end
%==========================================================================