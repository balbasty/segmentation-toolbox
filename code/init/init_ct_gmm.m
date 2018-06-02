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
    end
end

if ~isempty(pars_ct) && pars.do_template
    % There are data that is CT -> perform CT init routine (below)
    %----------------------------------------------------------------------
    
    M1  = numel(pars_ct.dat); 
    val = 0;
    
    % Calculate histogram for each CT image and store the histograms
    c   = {};
    x   = {};
    cnt = 1;
    mn  = Inf;
    mx  = -Inf;
    for m=1:M1 % Loop over CT populations
        S = pars_ct.dat{m}.S;
        if verbose, fprintf(1,'Getting histogram from subject (m=%d, S=%d) -    ',m,S); end
        for s=1:S % Loop over subjects
            
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
        
            % Get and mask image data
            Nii = nifti(pars_ct.dat{m}.V{s}.fname);
            img = single(Nii.dat(:,:,:));
            msk = spm_misc('msk_modality',img,'CT');
            img = img(msk(:));

            % Get histogram bins and create histogram
            x1 = min(img(:)):max(img(:));
            c1 = hist(img(:),x1);

            % Store global min and max values
            if max(x1)>mx, mx = max(x1); end
            if min(x1)<mn, mn = min(x1); end

            % Store counts and values
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
    for i=1:numel(c)
        x1 = x{i};        
        x1 = x1 + (abs(mn) + 1);
        
        msk     = false(1,numel(X));
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
%     nmn = -1040;
%     nmx = 3000;
%     nix = find(X>=nmn & X<=nmx);
%     X   = X(nix);
%     C   = C(nix);

    % For debugging
    SHOW_FIT = false;        
    
    % Fit VB-GMM to background
    k1              = 3;
    [~,ix]          = find(X<(-200 + val));
    [~,m1,b1,n1,W1] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k1,true,1e-8,false,SHOW_FIT);      
    lkp1            = ones(1,k1);
    
    % Fit VB-GMM to soft tissue
    k2              = 1;
    [~,ix]          = find(X>=(-200 + val) & X<(0 + val));
    [~,m2,b2,n2,W2] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k2,true,1e-8,false,SHOW_FIT);      
    lkp2            = 2*ones(1,k2);
    
    % For brain tissue, the parameters are set hard-coded 
    k3 = 7;
    m3 = [10 25 35 40 50 65 75] + val;
    b3 = ones(1,k3);
    n3 = ones(1,k3);
    W3 = ones(1,k3);
    lkp3 = [3 4 5 6 7 8 8];
    
    % Fit VB-GMM to bone
    k4              = 3;
    [~,ix]          = find(X>=(200 + val));
    [~,m4,b4,n4,W4] = spm_imbasics('fit_vbgmm2hist',C(ix),X(ix),k4,true,1e-8,false,SHOW_FIT);     
    lkp4            = Kb*ones(1,k4);
     
    m    = [m1,m2,m3,m4];
    b    = [b1,b2,b3,b4];
    n    = [n1,n2,n3,n4];
    W    = [W1,W2,W3,W4];
    part = [lkp1,lkp2,lkp3,lkp4];        
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
        
        gmm         = ngmm;
        gmm.pr.part = part; % just so that the partition can be loaded from the prior object
        part        = npart;
        mg          = nmg;
    end

    
    % Set CT GMM structs
    %----------------------------------------------------------------------    
    for m=1:M
        modality = pars.dat{m}.modality;
        if strcmp(modality,'CT')
            S = numel(pars.dat{m});
            for s=1:S
                pars.dat{m}.segment.mg       = mg;
                pars.dat{m}.segment.lkp.part = part;                
                pars.dat{m}.segment.gmm      = gmm;       
            end
        end
    end
end
%==========================================================================
