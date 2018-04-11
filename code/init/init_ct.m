function pars = init_ct(pars)
% Init GMM when data is CT
% FORMAT pars = init_ct(pars)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if ~isfield(pars,'K')
    pars.K = 22;
end

K = pars.K;

pars1 = read_images(pars,false);
M     = numel(pars1.dat);

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
        cnt             = cnt + 1;
    end
end

if ~isempty(npars1)
    h   = {};
    x   = {};
    cnt = 1;
    mn  = Inf;
    mx  = -Inf;
    for m=1:M
        S = numel(npars1.dat{m}.V);
        for s=1:S
            Nii = nifti(npars1.dat{m}.V{s}.fname);
            img = single(Nii.dat(:,:,:));
            msk = msk_modality(img,'CT');
            img = img(msk(:));

            x1     = min(img(:)):max(img(:));
            h1     = hist(img(:),x1);

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

    nmn   = -1050;
    ix_mn = find(X==nmn);
    X     = X(ix_mn:end);
    H     = H(ix_mn:end);

    [mg1,mn,vr1] = spm_imbasics('fit_gmm2hist',H,X,K);

    % Sort parameters according to mean
    [~,ix] = sort(mn);
    mn     = mn(ix);
%     mg1    = mg1(ix);
%     vr1    = vr1(ix);

    % Partition tissue classes
    p1 = mn<-200;
    p2 = mn>= -200 & mn<200;
    p3 = mn>=200;

    part = [ones(1,nnz(p1)) repelem(1:ceil(nnz(p2)/2),1,2) + 1 (2 + ceil(nnz(p2)/2))*ones(1,nnz(p3))];
    if numel(part)>K
        ix          = find(p2>0);
        part(ix(1)) = [];
    end
    Kb = max(part);

    % Reshape
    mn  = reshape(mn,1,K);
%     mg1 = reshape(mg1,1,K);
%     vr1 = reshape(vr1,1,1,K);
% 
%     sd = (mx - nmn)/K;
%     vr = sd^2;

    % Init GMM struct
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

    % Set pars struct
    pars.K = Kb;
    for m=1:M
        S = numel(npars1.dat{m});
        for s=1:S
            pars.dat{m}.segment.mg       = mg;
            pars.dat{m}.segment.lkp.part = part;
            pars.dat{m}.segment.gmm      = gmm;
        end
    end
end
%==========================================================================