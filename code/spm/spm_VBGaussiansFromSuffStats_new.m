function [out_L,fmom,mog] = spm_VBGaussiansFromSuffStats_new(mom,mog)

K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);
vb = mog.vb;

if vb
    % Variational Bayes solution
    
    nu0   = mog.pr.n;
    W0    = mog.pr.W;
    beta0 = mog.pr.b;
    m0    = mog.pr.m;
    
    if ~isempty(mog.po)
        nu   = mog.po.n;
        W    = mog.po.W;
        beta = mog.po.b;
        m    = mog.po.m;    
    else
        nu   = zeros(size(nu0));
        W    = zeros(size(W0));
        beta = zeros(size(beta0));
        m    = zeros(size(m0));
    end
    
    out_L = zeros(4,K);

    if N >= min(nu0)+1
        warning('SPM:Wishart','Can not normalise Wishart distribution (M=%d, nu_0=%f)', N,min(nu0));
    end        
else
    % Maximum likelihood solution
   
    vr0 = mog.vr0;
    
    if ~isempty(mog.m)
        mu = mog.m;
    else
        mu = zeros(N,K);
    end
    
    if ~isempty(mog.C)
        Sigma = mog.C;
    else
        Sigma = repmat(eye(N,N),[1 1 K]);        
    end
    
    C1 = 0;
    
    out_L = zeros(1,K);
end

fmom.s0 = zeros(1,K);
fmom.s1 = zeros(N,K);
fmom.S2 = zeros(N,N,K);

for k=1:K
    % Compute zeroeth moment
    s0 = 0;
    for i=1:numel(mom), s0 = s0 + mom(i).s0(k); end

    if vb && nargout==3
        nu(k)    = nu0(k)   + s0;
        beta(k)  = beta0(k) + s0;
    end

    L = -Inf;
    for iter=1:1024

        if vb
            muk    = m(:,k);
            Lambda = W(:,:,k)*nu(k);            
        else
            muk    = mu(:,k);
            Lambda = inv(Sigma(:,:,k));
        end

        s1     = zeros(N,1); % First moment
        S2     = zeros(N,N); % Second moment
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_L  = L;

        % E[log q(h)]
        L3 = 0;
        for i=1:numel(mom)
            if mom(i).s0(k)
                ind            = mom(i).ind;
                mux            = muk( ind,:);
                muy            = muk(~ind,:);
                Lambdayy       = Lambda(~ind,~ind);
                Lambdayx       = Lambda(~ind, ind);
                R              = Lambdayy\Lambdayx;
                Ex             = mom(i).s1(:,k)/mom(i).s0(k);
                Exx            = mom(i).S2(:,:,k)/mom(i).s0(k);
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = mom(i).s1(:,k);
                S2i( ind, ind) = mom(i).S2(:,:,k);
                s1i(~ind)      = mom(i).s0(k)*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = mom(i).s0(k)*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Lambdayy));
                S2i( ind,~ind) = mom(i).s0(k)*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;
                L3             = L3 + 0.5*mom(i).s0(k)*(logdet(Lambdayy) - size(Lambdayy,1)*(1+log(2*pi)));
            end
        end

        if vb
            if nargout==3
                m(:,k)   = (beta0(k)*m0(:,k) + s1)/beta(k);
                W(:,:,k) = inv(beta0(k)*m0(:,k)*m0(:,k)' + S2 - beta(k)*m(:,k)*m(:,k)' + inv(W0(:,:,k)));
                W(:,:,k) = threshold_eig(W(:,:,k));
            end
            
            % E[log p(g,h|mu,Lambda)]
            Eld = Elogdet(W(:,:,k),nu(k));
            L1  = 0.5*s0*(Eld - N*log(2*pi) - trace((m(:,k)*m(:,k)' + inv(beta(k)*nu(k)*W(:,:,k)) + S2/s0 - 2*m(:,k)*s1'/s0)*nu(k)*W(:,:,k)));

            % E[log p(mu,Lambda)]
            L2 = 0.5*N*(log(beta0/(2*pi)) - beta0/beta) + 0.5*(nu0-N)*Eld - logWishartDen(W0(:,:,k),nu0(k)) ...
               - 0.5*nu(k)*beta0(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)) - 0.5*nu(k)*trace(W0(:,:,k)\W(:,:,k));

            % E[log q(mu,Lambda)]
            L4 = 0.5*N*(log(beta(k)/(2*pi)) - 1 - nu(k)) + 0.5*(nu(k)-N)*Eld - logWishartDen(W(:,:,k),nu(k));

            % Negative Free Energy
            L = L1 + L2 - L3 - L4;
           %fprintf('%d\t%g\t%g\t%g\t%g\t%g\n', iter,L1, L2, L3, L4, L);                        
        else
            if nargout==3
                mu(:,k)      = s1/s0;
                Sigma(:,:,k) = (S2 - mu(:,k)*mu(:,k)'/s0 + N*vr0)/(s0 + N);  
            end
            
            L = -0.5*s0*(logdet(Sigma(:,:,k)) + N*log(2*pi)) ...
                -0.5*s0*trace((mu(:,k)*mu(:,k)' + S2/s0 - 2*mu(:,k)*s1'/s0)/Sigma(:,:,k)) ...
                -L3;
           %fprintf('%d\t%g\n', iter,L);
        end
        if L-old_L < s0*N*1e-10, break; end
    end
    
    fmom.s0(k)     = s0;
    fmom.s1(:,k)   = s1;
    fmom.S2(:,:,k) = S2;
    
    if vb
        out_L(:,k) = [L1 L2 L3 L4]';        
    else
        C1 = C1 + (S2 - mu(:,k)*mu(:,k)'/s0);
        
        out_L(k)   = L;
    end
end

out_L = sum(sum(out_L));

mog.s0 = fmom.s0;
if vb
    % Posterior distributions
    po.m   = m;
    po.W   = W;
    po.n   = nu;
    po.b   = beta;
    mog.po = po;    
else
    % Maximum likelihood estimates   
    vr1 = zeros(N,N,K);
    for k=1:K
        vr1(:,:,k) = (C1 + N*vr0)/(sum(fmom.s0) + N);
    end 
    
    mog.vr1 = vr1;
    mog.m   = mu;
    mog.C   = Sigma;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = Elogdet(W,nu)
% E[log(det(Lambda))], where Lambda ~ W(W,nu)
M = size(W,1);
e = logdet(W) + M*log(2);
for m=1:M
    e = e+psi((nu-m+1)/2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = logdet(W)
% log(det(W))
l = 2*sum(log(diag(chol(W))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lZ = logWishartDen(W,nu)
% Log of normalising term of log W(W,nu)
M     = size(W,1);
if M >= nu+1
   %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
    lZ = 0;
    return;
end
lGamM = M*(M-1)/4*log(pi);
for m=1:M
    lGamM = lGamM+gammaln((nu+1-m)/2);
end
lZ    = 0.5*nu*(logdet(W) + M*log(2)) + lGamM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%