function pars = init_ct(pars,test_level)
% Init GMM when data is CT
% FORMAT pars = init_ct(pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, test_level = 0; end

if ~isfield(pars,'K')
    pars.K = 10;
end
Kb = pars.K; % Requested number of tissue classes

% Read only CT data-sets into a new pars structure (pars1)
pars1  = read_images(pars,false);
M      = numel(pars1.dat);

if test_level
    % Adjust number of subjects if testing
    for m=1:M
        if     test_level==2 || test_level==3, pars1.dat{m}.S = min(8,pars1.dat{m}.S);
        elseif test_level==1,                  pars1.dat{m}.S = 1;   
        end 
    end
end

tot_S  = 0;
npars1 = {};
cnt    = 1;
for m=1:M
    if isfield(pars.dat{m},'modality')
        modality = pars.dat{m}.modality;
    else
        continue
    end
    if strcmp(modality,'CT')
        npars1.dat{cnt} = pars1.dat{m};
        tot_S           = npars1.dat{cnt}.S;
        cnt             = cnt + 1;
    end
end
clear pars1

if ~isempty(npars1) && tot_S>1
    % There are data that is CT -> perform CT init routine (below)
    %----------------------------------------------------------------------
    
    % Calculate histogram for each CT image and store the histograms
    h   = {};
    x   = {};
    cnt = 1;
    mn  = Inf;
    mx  = -Inf;
    for m=1:M
        S = npars1.dat{m}.S;
        for s=1:S
            Nii = nifti(npars1.dat{m}.V{s}.fname);
            img = single(Nii.dat(:,:,:));
            msk = msk_modality(img,'CT');
            img = img(msk(:));

            x1 = min(img(:)):max(img(:));
            h1 = hist(img(:),x1);

            if max(x1)>mx
                mx = max(x1);
            end

            if min(x1)<mn
                mn = min(x1);
            end

            h{cnt} = single(h1);
            x{cnt} = single(x1);

            cnt = cnt + 1;
        end
    end
    clear h1 x1

    % Sum all histograms into (H)
    X   = mn:mx;
    H   = zeros(1,numel(X));
    msk = false(1,numel(X));
    for i=1:numel(h)
        x1      = x{i};
        x1      = x1 + abs(mn) + 1;
        msk(x1) = true;

        H(msk) = H(msk) + h{i};

        msk(1,:) = false;
    end
    clear h x

    % Remove intensity values smaller than a threshold (nmn)
    nmn   = -1050;
    ix_mn = find(X>=nmn);
    X     = X(ix_mn(1):end);
    H     = H(ix_mn(1):end);

    % Keep fitting a K1-component GMM to the accumulated histogram until a
    % Kb-component GMM can be constructed. Kb are the number of requested
    % tissue classes in the TPM to be built. For CT each tissue class are
    % represented by a specific number of Gaussians (see [1] below). 
    %----------------------------------------------------------------------
    for K1=(2*Kb + 2):50
        % Fit GMM to histogram
        [~,mn,vr] = spm_imbasics('fit_gmm2hist',H,X,K1);

        % Sort parameters according to mean
        [~,ix] = sort(mn);
        mn     = mn(ix);
        vr     = vr(ix);
        
        % 1. Partition tissue classes
        p1 = mn<-50; % Background intensities have nnz(p1) Gaussians
        p2 = mn>= -50 & mn<200; % Brain intensities have 2 Gaussians each
        p3 = mn>=200; % Bone intensities have nnz(p2) Gaussians

        % Create lkp.part partitioning
        part = [ones(1,nnz(p1)) repelem(1:ceil(nnz(p2)/2),1,2) + 1 (2 + ceil(nnz(p2)/2))*ones(1,nnz(p3))];
                
        if numel(part)>K1
            % It can happen that there are more elements in part than in K1,
            % the number of means then needs to be adjusted accordingly
            ix  = find(p2>0);
            vr0 = vr(ix(1));
            mn0 = mn(ix(1));
            w   = 1./(1 + exp(-0.25)) - 0.5;            
            mn1 = sqrtm(vr0)*randn(1,2)*w + repmat(mn0,[1,2]);
            
            mn = [mn(p1); mn1(1); mn1(2); mn(ix(2):end)];            
        end        
        
        if max(part)==Kb
            % The requested number of tissue classes have been obtained,
            % stop
            K = numel(mn);
            break 
        end
    end
    
    % Init CT GMM structs
    %----------------------------------------------------------------------
    mn  = reshape(mn,1,K);
    gmm = struct;
    mg  = zeros(1,K);
    for m=1:M
        S = numel(npars1.dat{m});
        for s=1:S
            for k=1:Kb
                ix = find(part==k);

                mg(ix) = 1/nnz(ix);

                gmm.pr.m(:,ix)   = mn(:,ix);
                gmm.pr.b(ix)     = 1;
                gmm.pr.n(ix)     = 1;
                gmm.pr.W(:,:,ix) = 1;

                gmm.po.m(:,ix)   = mn(:,ix);
                gmm.po.b(ix)     = 1;
                gmm.po.n(ix)     = 1;
                gmm.po.W(:,:,ix) = 1;                
            end
        end
    end

    % Set CT GMM structs
    %----------------------------------------------------------------------
    pars.K = Kb;
    for m=1:M
        modality = pars.dat{m}.modality;
        if strcmp(modality,'CT')
            S = numel(npars1.dat{m});
            for s=1:S
                pars.dat{m}.segment.mg       = mg;
                pars.dat{m}.segment.lkp.part = part;
                pars.dat{m}.segment.gmm      = gmm;
            end
        end
    end
end
%==========================================================================