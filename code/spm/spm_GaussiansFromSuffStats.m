function [lb,lh,mog,s0] = spm_GaussiansFromSuffStats(mom,mog)

tiny = eps*eps;
K    = size(mom(1).s0,2);
N    = numel(mom(1).ind);
vb   = mog.vb;

if vb
    % Variational Bayes solution
    %----------------------------------------------------------------------
    
    if isfield(mog,'pr') && isfield(mog,'po')
        n0 = mog.pr.n;
        W0 = mog.pr.W;
        b0 = mog.pr.b;
        m0 = mog.pr.m;
        if numel(n0)  == 1, n0 = repmat(n0,[1 K  ]); end
        if size(W0,3) == 1, W0 = repmat(W0,[1 1 K]); end
        if numel(b0)  == 1, b0 = repmat(b0,[1 K  ]); end
        if size(m0,2) == 1, m0 = repmat(m0,[1 K  ]); end
    
        n = mog.po.n;
        W = mog.po.W;
        b = mog.po.b;
        m = mog.po.m;  
        
        vr1 = [];
    elseif isfield(mog,'pr') && ~isfield(mog,'po')
        n0 = mog.pr.n;
        W0 = mog.pr.W;
        b0 = mog.pr.b;
        m0 = mog.pr.m;
        if numel(n0)  == 1, n0 = repmat(n0,[1 K  ]); end
        if size(W0,3) == 1, W0 = repmat(W0,[1 1 K]); end
        if numel(b0)  == 1, b0 = repmat(b0,[1 K  ]); end
        if size(m0,2) == 1, m0 = repmat(m0,[1 K  ]); end
        
        n  = n0;
        W  = W0;
        b  = b0;
        m  = m0;
        
        vr1 = [];
    else
        m0 = zeros(N,K);
        b0 = ones(1,K);
        n0 = (N - .999)*ones(1,K);
        W0 = repmat(eye(N,N),[1 1 K]);
        
        n  = n0;
        W  = W0;
        b  = b0;
        m  = m0;
        
        vr1 = eye(N,N);
    end
    
    if N >= min(n0)+1
        warning('SPM:Wishart','Can not normalise Wishart distribution (M=%d, nu_0=%f)', N,min(n0));
    end        
else
    % Maximum likelihood solution       
    %----------------------------------------------------------------------
    
    if isfield(mog,'mn') && isfield(mog,'vr')
        mn  = mog.mn;   
        vr  = mog.vr; 
        
        vr1 = [];
    else
        mn  = zeros(N,K);
        vr  = repmat(eye(N,N),[1 1 K]);
        
        vr1 = eye(N,N);
    end
end

s0 = zeros(1,K);
s1 = zeros(N,K);
S2 = zeros(N,N,K);
lb = 0;
lh = 0;
for k=1:K
    % Compute zeroeth moment
    s0(k) = 0;
    for i=1:numel(mom), s0(k) = s0(k) + mom(i).s0(k); end

    if vb && nargout>=3
        n(k) = n0(k) + s0(k);
        b(k) = b0(k) + s0(k);
    end

    Lk = -Inf;
    for iter=1:1024
        oLk  = Lk;

        if vb
            mnk     = m(:,k);
            Lambdak = W(:,:,k)*n(k);            
        else
            mnk     = mn(:,k);
            Lambdak = inv(vr(:,:,k));
        end

        s1(:,k)   = zeros(N,1); % First moment
        S2(:,:,k) = zeros(N,N); % Second moment
        s1i       = zeros(N,1);
        S2i       = zeros(N,N);
        
        % E[log q(h)]
        L3 = 0;
        for i=1:numel(mom)
            if mom(i).s0(k)>tiny
                ind            = mom(i).ind;
                mnkx           = mnk( ind,:);
                mnky           = mnk(~ind,:);
                Lambdayy       = Lambdak(~ind,~ind);
                Lambdayx       = Lambdak(~ind, ind);
                R              = Lambdayy\Lambdayx;
                Ex             = mom(i).s1(:,k)/mom(i).s0(k);
                Exx            = mom(i).S2(:,:,k)/mom(i).s0(k);
                tmp            = R*(mnkx-Ex)*mnky';
                s1i( ind)      = mom(i).s1(:,k);
                S2i( ind, ind) = mom(i).S2(:,:,k);
                s1i(~ind)      = mom(i).s0(k)*(R*(mnkx-Ex) + mnky);
                S2i(~ind,~ind) = mom(i).s0(k)*(R*(mnkx*mnkx'-mnkx*Ex'-Ex*mnkx'+Exx)*R' + tmp + tmp' + mnky*mnky' + inv(Lambdayy));
                S2i( ind,~ind) = mom(i).s0(k)*(R*(mnkx*Ex'-Exx) + mnky*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1(:,k)        = s1(:,k) + s1i;
                S2(:,:,k)      = S2(:,:,k) + S2i;
                L3             = L3 + 0.5*mom(i).s0(k)*(LogDet(Lambdayy) - size(Lambdayy,1)*(1+log(2*pi)));
            end
        end
        
        if vb            
            mom1.s0 = s0(k);
            mom1.s1 = s1(:,k);
            mom1.S2 = S2(:,:,k);                
            mom1    = mom_John2Bishop(mom1);  
                
            if nargout>=3
                m(:,k) = (b0(k)*m0(:,k) + s0(k)*mom1.s1)./b(k);
                
                mlt1     = b0(k).*s0(k)/(b0(k) + s0(k));
                diff1    = mom1.s1 - m0(:,k);
                W(:,:,k) = inv(inv(W0(:,:,k)) + s0(k)*mom1.S2 + mlt1*(diff1*diff1'));                
                W(:,:,k) = threshold_eig(W(:,:,k));
            end
            
            % E[log p(g,h|mu,Lambda)]
            Eld = Elogdet(W(:,:,k),n(k));
            L1  = 0.5*s0(k)*(Eld - N*log(2*pi) - trace((m(:,k)*m(:,k)' + inv(b(k)*n(k)*W(:,:,k)) + S2(:,:,k)/s0(k) - 2*m(:,k)*s1(:,k)'/s0(k))*n(k)*W(:,:,k)));

            % E[log p(mu,Lambda)]
            L2 = 0.5*N*(log(b0(k)/(2*pi)) - b0(k)/b(k)) + 0.5*(n0(k)-N)*Eld - logWishartDen(W0(:,:,k),n0(k)) ...
                 - 0.5*n(k)*b0(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)) - 0.5*n(k)*trace(W0(:,:,k)\W(:,:,k));

            % E[log q(mu,Lambda)]
            L4 = 0.5*N*(log(b(k)/(2*pi)) - 1 - n(k)) + 0.5*(n(k)-N)*Eld - logWishartDen(W(:,:,k),n(k));

            % Negative Free Energy
            Lk = L1 + L2 - L3 - L4;
           %fprintf('%d\t%g\t%g\t%g\t%g\t%g\n', iter,L1, L2, L3, L4, L);                        
        else
            if nargout>=3
                mn(:,k) = s1(:,k)/s0(k);
                vr(:,:,k) = (S2(:,:,k) - s1(:,k)*s1(:,k)'/s0(k) + N*mog.vr0)/(s0(k) + N);         
            end
            
            Lk = - 0.5*s0(k)*(LogDet(vr(:,:,k)) + N*log(2*pi)) ...
                 - 0.5*s0(k)*trace((mn(:,k)*mn(:,k)' + S2(:,:,k)/s0(k) - 2*mn(:,k)*s1(:,k)'/s0(k))/vr(:,:,k)) ...
                 - L3;
%            fprintf('%d\t%g\n', iter,L);
        end
        if Lk-oLk < s0(k)*N*1e-10
            break; 
        end
    end
       
    if ~isempty(vr1)
        vr1 = vr1 + (S2(:,:,k) - s1(:,k)*s1(:,k)'/s0(k));   
    end
    
    if vb && nargout<3  
        po.m = m(:,k);
        po.b = b(k);
        po.W = W(:,:,k);
        po.n = n(k);

        pr.m = m0(:,k);
        pr.b = b0(k);
        pr.W = W0(:,:,k);
        pr.n = n0(k);

        lh = lh - L3;
        lb = lb + lowerbound_GaussWishart(mom1,po,pr);
    end
end

if vb
    if isfield(mog,'pr') && isfield(mog,'po')
        mog.po.m = m;
        mog.po.W = W;
        mog.po.n = n;
        mog.po.b = b;
    else
        if ~isempty(vr1)
            vr1 = (vr1 + N*mog.vr0)/(sum(s0) + N);
        end
        
        mom    = [];
        mom.s0 = s0;
        mom.s1 = s1;
        mom.S2 = S2;
        mom    = mom_John2Bishop(mom);  
        
        if ~isfield(mog,'pr')
            mog.pr.m = zeros(N,K);
            mog.pr.b = ones(1,K);
            mog.pr.n = (N - .999)*ones(1,K);                                   
            mog.pr.W = repmat(eye(N,N),1,1,K);
%             mog.pr.m = mom.s1;
%             mog.pr.b = 1e-3*s0;
%             mog.pr.n = 1e-3*s0 + N - 1;                                   
%             mog.pr.W = zeros(N,N,K);
%             for k=1:K
%                 mog.pr.W(:,:,k) = inv(vr1)/mog.pr.n(k);
%             end
        end
        
        mog.po.n = mog.pr.n + s0;
        mog.po.b = mog.pr.b + s0;            
        for k=1:K
            mog.po.m(:,k)   = (mog.pr.b(k)*mog.pr.m(:,k) + s0(k)*mom.s1(:,k))./mog.po.b(k);
            mlt1            = mog.pr.b(k).*s0(k)/(mog.pr.b(k) + s0(k));
            diff1           = mom.s1(:,k) - mog.pr.m(:,k);
            mog.po.W(:,:,k) = inv(inv(mog.pr.W(:,:,k)) + s0(k)*mom.S2(:,:,k) + mlt1*(diff1*diff1')); 
            mog.po.W(:,:,k) = threshold_eig(mog.po.W(:,:,k));
        end
    end
else
    mog.mn  = mn;
    if isempty(vr1)
        mog.vr  = vr;
    else
        vr1 = (vr1 + N*mog.vr0)/(sum(s0) + N);
        for k=1:K
            mog.vr(:,:,k) = vr1;
        end
    end
end
%==========================================================================