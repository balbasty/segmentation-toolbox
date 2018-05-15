function [out_L,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm)
% Gaussians from sufficient statistics (with missing data)
%
% FORMAT [post,L] = spm_VBGaussiansFromSuffStats(mom,priors)
% mom        - Sufficient statistics from spm_SuffStats
% priors     - A structure having fields
%              nu   (1xK)
%              W    (MxMxK)
%              beta (1xK)
%              m    (MxK)
%              These definine priors for K Normal-Wishart
%              distributions:
%              N(mu|m,inv(beta Lambda)) W(Lambda|W,nu)
% post       - A structure having the same fields as the priors,
%              which define the posterior distributions for
%              K Normal-Wishart distributions.
% L            A 4xK matrix, encoding the expecatations used
%              for computing Free Energy (where g are observed
%              data, and h are missing).
%              * L(:,1) = E[log p(g,h|mu,Lambda)]
%              * L(:,2) = E[log p(mu,Lambda)]
%              * L(:,3) = E[log q(h)]
%              * L(:,4) = E[log q(mu,Lambda)]
%
% This is the Bayesian version.  See Bishop's ``Pattern Recognition
% and Machine Learning'' textbook for further elaborations.
% The overall Negative Free Energy objective function is
% sum(L(:,1) + L(:,2) - L(:,3) - L(:,4))
%
%
% FORMAT [N,mu,Sigma,L] = spm_VBGaussiansFromSuffStats(mom)
% mom        - Sufficient statistics from spm_SuffStats
% N          - number of samples
% mu         - maximum likelihood estimates of mu (MxK)
% Sigma      - maximum likelihood estimates of inv(Lambda) (MxMxK)
% L          - log-likelihood
%
% This is the maximum likelihood version, which aims to outperform
% MATLAB's nanmean(X) and nancov(X,1).
%
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

K    = size(mom(1).s0,2);
N    = numel(mom(1).ind);
tiny = eps*eps;

nu0   = gmm.pr.n;
W0    = gmm.pr.W;
if isfield(gmm.pr, 'ldW')
    ldW0  = gmm.pr.ldW;
else
    ldW0 = zeros(1,K);
    for k=1:K
        ldW0(k) = spm_matcomp('LogDet', W0(:,:,k));
    end
end
beta0 = gmm.pr.b;
m0    = gmm.pr.m;
if numel(nu0)   == 1, nu0   = repmat(  nu0,[1 K  ]); end
if size(W0,3)   == 1, W0    = repmat(   W0,[1 1 K]); end
if numel(beta0) == 1, beta0 = repmat(beta0,[1 K  ]); end
if size(m0,2)   == 1, m0    = repmat(   m0,[1 K  ]); end

% Init posteriors
%----------------------------------------------------------------------
if isfield(gmm,'po') 
    nu    = gmm.po.n;
    W     = gmm.po.W;
    beta  = gmm.po.b;
    m     = gmm.po.m;
else                  
    m    = gmm.pr.m;                       
    beta = gmm.pr.b;
    nu   = gmm.pr.n; 
    W    = gmm.pr.W;        
end

out_L = zeros(4,K);

if N >= min(nu0)+1
    warning('SPM:Wishart','Can not normalise Wishart distribution (M=%d, nu_0=%f)', N,min(nu0));
end

for k=1:K
    % Compute zeroeth moment
    s0 = 0;
    for i=2:numel(mom), s0 = s0 + mom(i).s0(k); end

    if nargout>1
        nu(k)   = nu0(k)   + s0;
        beta(k) = beta0(k) + s0;
    end

    L = -Inf;
    for iter=1:1024

        muk    = m(:,k);
        Lambda = W(:,:,k)*nu(k);

        s1     = zeros(N,1); % First moment
        S2     = zeros(N,N); % Second moment
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_L  = L;

        % E[log q(h)]
        L3     = 0;
        for i=1:numel(mom)
            if mom(i).s0(k)>tiny
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
                L3             = L3 + 0.5*mom(i).s0(k)*(spm_matcomp('LogDet',Lambdayy) - size(Lambdayy,1)*(1+log(2*pi)));
            end
        end

        if nargout>1
            m(:,k)   = (beta0(k)*m0(:,k) + s1)/beta(k);

            Wk = spm_matcomp('Inv',beta0(k)*m0(:,k)*m0(:,k)' + S2 - beta(k)*m(:,k)*m(:,k)' + spm_matcomp('Inv',W0(:,:,k)));
            Wk = spm_matcomp('threshold_eig',Wk);
            [~,p] = chol(Wk);
            if p == 0
                W(:,:,k) = Wk;
            end
        end

        % E[log p(g,h|mu,Lambda)]
        Eld = spm_prob('Wishart','Elogdet',W(:,:,k),nu(k));
        L1  = 0.5*s0*(Eld - N*log(2*pi) - trace((m(:,k)*m(:,k)' + spm_matcomp('Inv',beta(k)*nu(k)*W(:,:,k)) + S2/s0 - 2*m(:,k)*s1'/s0)*nu(k)*W(:,:,k)));

        % E[log p(mu,Lambda)]
        L2 = 0.5*N*(log(beta0(k)/(2*pi)) - beta0(k)/beta(k)) + 0.5*(nu0(k) - N)*Eld ...
            - nu0(k)*N/2*log(2) - nu0(k)/2*ldW0(k) - spm_prob('LogGamma',nu0(k)/2,N) ...
            - 0.5*nu(k)*beta0(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)) - 0.5*nu(k)*trace(W0(:,:,k)\W(:,:,k));
            
        % E[log q(mu,Lambda)]
        L4 = 0.5*N*(log(beta(k)/(2*pi)) - 1 - nu(k)) + 0.5*(nu(k)-N)*Eld ...
            - spm_prob('Wishart','logZ',W(:,:,k),nu(k));

        % Negative Free Energy
        L = L1 + L2 - L3 - L4;
       %fprintf('%d\t%g\t%g\t%g\t%g\t%g\n', iter,L1, L2, L3, L4, L);
        
        if L-old_L < s0*N*1e-10, break; end
    end    
        
    out_L(:,k) = [L1 L2 L3 L4]';
end

% Posterior distributions
gmm.po.m = m;
gmm.po.W = W;
gmm.po.n = nu;
gmm.po.b = beta;
%==========================================================================