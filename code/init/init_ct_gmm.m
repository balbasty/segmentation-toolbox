function pars = init_ct_gmm(pars,verbose)
% Init GMM when data is CT
% FORMAT pars = init_ct_gmm(pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin<2, verbose = true; end

Kb = pars.K; 
M  = numel(pars.dat);

% Read CT populations into new struct (pars_ct)
S_ct    = 0;
pars_ct = {};
cnt     = 1;
healthy = []; 
for m=1:M
    if isfield(pars.dat{m},'modality')
        modality = pars.dat{m}.modality;
    else
        continue
    end
    if strcmp(modality,'CT')
        pars_ct.dat{cnt} = pars.dat{m};
        S_ct             = S_ct + pars_ct.dat{cnt}.S;
        cnt              = cnt + 1;
        
        healthy = [healthy pars.dat{m}.healthy]; 
    end
end

if ~isempty(pars_ct) && pars.do_template
    % There are data that is CT -> perform CT init routine (below)
    %----------------------------------------------------------------------
    
    M1            = numel(pars_ct.dat); 
    ngauss_lesion = pars_ct.dat{1}.segment.ct.ngauss_lesion; 
    ngauss_brain  = pars_ct.dat{1}.segment.ct.ngauss_brain;     
    
    if sum(healthy)==M1 
        % All healthy 
        flag_healthy_ct = 0; 
    elseif sum(healthy)==0 
        % All lesion 
        flag_healthy_ct = 1; 
    else 
        % Mix of healthy and lesion         
        flag_healthy_ct = 2; 
    end 
    
    % Calculate histogram for each CT image and store the histograms
    c   = {};
    x   = {};
    cnt = 1;
    mn  = Inf;
    mx  = -Inf;
    for m=1:M1
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
    C = C/S_ct;
        
    % Smooth histogram a little bit
%     krn = spm_smoothkern(4,(-256:256)',0);
%     C   = conv(C,krn,'same');
    
    % Remove intensity values smaller than an upper and lower threshold
    nmn = -1040;
    nmx = 3000;
    nix = find(X>=nmn & X<=nmx);
    X   = X(nix);
    C   = C(nix);

    % k1 = 5, k2 = 1, ngauss_brain = 1 | old was k1 = 4, k2 = 2,
    % ngauss_brain = 2
    
    % Fit VB-GMM to background
    k1              = 5;
    [~,ix]          = find(X<-200);
    [~,m1,b1,n1,W1] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k1);      
    lkp1            = ones(1,k1);
    
    % Fit VB-GMM to soft tissue
    k2              = 1;
    [~,ix]          = find(X>=-200 & X<0);
    [~,m2,b2,n2,W2] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k2);      
    lkp2            = 2*ones(1,k2);
    
    % Fit VB-GMM to brain    
    if flag_healthy_ct==0 
        k3 = 5;
    elseif flag_healthy_ct==1 || flag_healthy_ct==2 
        k3 = 6;
    end
    [~,ix]          = find(X>=0 & X<80);
    [~,m3,b3,n3,W3] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k3);     
    lkp3            = repelem(1:k3,1,ngauss_brain);
        
    if flag_healthy_ct~=0
        [~,ix_les] = min(abs(m3 - 60)); % 60 should in CT contain (mostly) lesion   
        
        lkp3(lkp3==ix_les) = []; 
        ix0                = find(lkp3==ix_les - 1); 
        ix1                = find(lkp3==ix_les + 1); 
        if isempty(ix0) 
            lkp3 = [ix_les*ones(1,ngauss_lesion) lkp3(ix1(1):end)];     
        elseif isempty(ix1) 
            lkp3 = [lkp3(1:ix0(end)) ix_les*ones(1,ngauss_lesion)]; 
        else 
            lkp3 = [lkp3(1:ix0(end)) ix_les*ones(1,ngauss_lesion) lkp3(ix1(1):end)]; 
        end     
        
        ix_les = ix_les + 2;
    end
    
    tmp_gmm.pr.m = m3;
    tmp_gmm.pr.n = n3;
    tmp_gmm.pr.b = b3;
    tmp_gmm.pr.W = W3;
    tmp_gmm.po   = tmp_gmm.pr;
    tmp_gmm      = more_gaussians(tmp_gmm,lkp3);
    m3           = tmp_gmm.po.m;
    n3           = tmp_gmm.po.n;
    b3           = tmp_gmm.po.b;
    W3           = squeeze(tmp_gmm.po.W)';
    clear tmp_gmm
    
    lkp3 = lkp3 + 2;
    
    % Fit VB-GMM to calcification        
    k4              = 1;
    [~,ix]          = find(X>=80 & X<200);
    [~,m4,b4,n4,W4] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k4); 
    lkp4            = (Kb - 1)*ones(1,k4);
    
    % Fit VB-GMM to bone
    k5              = 3;
    [~,ix]          = find(X>=200);
    [~,m5,b5,n5,W5] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k5);     
    lkp5            = Kb*ones(1,k5);
     
    m    = [m1,m2,m3,m4,m5];
    b    = [b1,b2,b3,b4,b5];
    n    = [n1,n2,n3,n4,n5];
    W    = [W1,W2,W3,W4,W5];
    part = [lkp1,lkp2,lkp3,lkp4,lkp5];        
    K    = numel(part);       
    
    % Init VB-GMM posteriors
    %----------------------------------------------------------------------    
    gmm = struct;
    mg  = zeros(1,K);
    for k=1:Kb
        ix = find(part==k);

        mg(ix) = 1/nnz(ix);

        gmm.pr.m(:,ix)   = m(ix);
        gmm.pr.b(ix)     = b(ix);
        gmm.pr.n(ix)     = n(ix);
        gmm.pr.W(:,:,ix) = reshape(W(ix),1,1,numel(ix));                    
    end
    
    % If user-specified class order is given, change order of GMM
    % parameters
    %----------------------------------------------------------------------    
    class_ix = pars_ct.dat{1}.segment.class_ix; 
    if ~isempty(class_ix)

        cnt   = 1;
        npart = zeros(1,K);
        nmg   = zeros(1,K);
        for k=1:Kb
            kk = find(class_ix(k)==part);
            for k1=kk
                ngmm.pr.b(cnt)     = gmm.pr.b(k1);
                ngmm.pr.n(cnt)     = gmm.pr.n(k1);
                ngmm.pr.m(:,cnt)   = gmm.pr.m(:,k1);
                ngmm.pr.W(:,:,cnt) = gmm.pr.W(:,:,k1);
                npart(cnt)         = k;
                nmg(cnt)           = 1/numel(kk);
                
                cnt = cnt + 1;                    
            end
        end
        
        gmm  = ngmm;
        part = npart;
        mg   = nmg;
    end

    
    % Set CT GMM structs
    %----------------------------------------------------------------------    
    for m=1:M
        modality = pars.dat{m}.modality;
        healthy  = pars.dat{m}.healthy;
        if strcmp(modality,'CT')
            S = numel(pars.dat{m});
            for s=1:S
                pars.dat{m}.segment.mg       = mg;
                pars.dat{m}.segment.lkp.part = part;                
                pars.dat{m}.segment.gmm      = gmm;
                
                if healthy && flag_healthy_ct==2
                    % If labelleled healthy, remove class with intensity
                    % closest to intensity of blood in CT.
                    pars.dat{m}.segment.lkp.rem  = ix_les;
                end                
            end
        end
    end
end
%==========================================================================
