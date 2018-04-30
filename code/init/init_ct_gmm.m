function pars = init_ct_gmm(pars,verbose)
% Init GMM when data is CT
% FORMAT pars = init_ct_gmm(pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin<2, verbose = true; end

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
        tot_S            = tot_S + pars_ct.dat{cnt}.S;
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
        if verbose, fprintf(1,'Getting histogram from subject (m=%d, S=%d) -    ',m,S); end
        for s=1:S
            
            if verbose                
                if s<10
                    fprintf(1,'\b%d',s);
                elseif s<100
                    fprintf(1,'\b\b%d',s);
                elseif s<1000
                    fprintf(1,'\b\b\b%d',s);
                elseif s<10000
                    fprintf(1,'\b\b\b\b%d',s);
                end                
            end
        
            Nii = nifti(pars_ct.dat{m}.V{s}.fname);
            img = single(Nii.dat(:,:,:));
            msk = spm_misc('msk_modality',img,'CT');
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
        fprintf(' | DONE!\n')   
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

    % Divide count by total number of subjects
    C = C/tot_S;
        
    % Smooth histogram a little bit
%     krn = spm_smoothkern(4,(-256:256)',0);
%     C   = conv(C,krn,'same');
    
    % Remove intensity values smaller than an upper and lower threshold
    nmn = -1040;
    nmx = 3000;
    nix = find(X>=nmn & X<=nmx);
    X   = X(nix);
    C   = C(nix);

    % Fit VB-GMM to background
    [~,ix]          = find(X<-200);
    [~,m1,b1,n1,W1] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),4);      
    lkp1            = ones(1,4);
    
    % Fit VB-GMM to soft tissue
    [~,ix]          = find(X>=-200 & X<0);
    [~,m2,b2,n2,W2] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),2);      
    lkp2            = 2*ones(1,2);
    
    % Fit VB-GMM to brain    
    [~,ix]          = find(X>=0 & X<80);
    [~,m3,b3,n3,W3] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),6); 
    lkp3            = 3:numel(m3) + 2;            
        
    % Fit VB-GMM to bone
    [~,ix]          = find(X>=80);
    [~,m4,b4,n4,W4] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),3);     
    lkp4            = Kb*ones(1,3);
     
    m    = [m1,m2,m3,m4];
    b    = [b1,b2,b3,b4];
    n    = [n1,n2,n3,n4];
    W    = [W1,W2,W3,W4];
    part = [lkp1,lkp2,lkp3,lkp4];      
  
    % Get index of lesions class
    [~,ix_rem] = min(abs(m - 60)); % 60 should in CT contain (mostly) lesion
    
    % Model lesion with a GMM
    ngauss_lesion = pars_ct.dat{1}.segment.ct.ngauss_lesion;
    vr_l          = 1./(n(ix_rem).*W(ix_rem));    
    w_l           = 1./(1 + exp(-(ngauss_lesion - 1)*0.25)) - 0.5;
    pr_l          = 1./(vr_l*(1 - w_l));  
    
    m_l = sqrtm(vr_l)*randn(1,ngauss_lesion)*w_l + repmat(m(ix_rem),[1,ngauss_lesion]);      
    b_l = b(ix_rem)/ngauss_lesion*ones(1,ngauss_lesion);    
    n_l = n(ix_rem)/ngauss_lesion*ones(1,ngauss_lesion);    
    W_l = pr_l/n(ix_rem)*ones(1,ngauss_lesion);                
    
    % Adjust GMM parameters
    m = [m(1:ix_rem - 1) m_l m(ix_rem + 1:end)];
    b = [b(1:ix_rem - 1) b_l b(ix_rem + 1:end)];
    n = [n(1:ix_rem - 1) n_l n(ix_rem + 1:end)];
    W = [W(1:ix_rem - 1) W_l W(ix_rem + 1:end)];        
    
    part = [part(1:ix_rem - 1) repmat(part(ix_rem),1,ngauss_lesion) part(ix_rem + 1:end)];
    K    = numel(part);   
    
    ix_rem = ix_rem:ix_rem + (ngauss_lesion - 1);
    
    % Init VB-GMM posteriors
    %----------------------------------------------------------------------    
    gmm = struct;
    mg  = zeros(1,K);
    for m1=1:M
        S = numel(pars_ct.dat{m1});
        for s=1:S
            for k=1:Kb
                ix = find(part==k);

                mg(ix) = 1/nnz(ix);

                gmm.pr.m(:,ix)   = m(ix);
                gmm.pr.b(ix)     = b(ix);
                gmm.pr.n(ix)     = n(ix);
                gmm.pr.W(:,:,ix) = reshape(W(ix),1,1,numel(ix));

                gmm.po.m(:,ix)   = m(ix);
                gmm.po.b(ix)     = b(ix);
                gmm.po.n(ix)     = n(ix);
                gmm.po.W(:,:,ix) = reshape(W(ix),1,1,numel(ix));      
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
                    pars.dat{m}.segment.lkp.rem  = part(ix_rem(1));
                    
                    pars.dat{m}.segment.gmm.po.m(ix_rem) = 0;
                    pars.dat{m}.segment.gmm.pr.m(ix_rem) = 0;
                    pars.dat{m}.segment.gmm.po.n(ix_rem) = 1;
                    pars.dat{m}.segment.gmm.pr.n(ix_rem) = 1;
                    pars.dat{m}.segment.gmm.po.b(ix_rem) = 1;
                    pars.dat{m}.segment.gmm.pr.b(ix_rem) = 1;
                    pars.dat{m}.segment.gmm.po.W(ix_rem) = 1;
                    pars.dat{m}.segment.gmm.pr.W(ix_rem) = 1;
                end
            end
        end
    end
end
%==========================================================================
