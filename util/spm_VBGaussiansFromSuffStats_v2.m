function varargout = spm_VBGaussiansFromSuffStats_v2(mom,priors)
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

K  = size(mom(1).s0,2);
M  = numel(mom(1).ind);

if nargin < 2,
    % Maximum likelihood solution
    ml    = true;
    N     = zeros(1,K);
    mu    = zeros(M,K);
    Sigma = repmat(eye(M,M),[1 1 K]);
    out_L = zeros(1,K);
else
    % Variational Bayes solution
    ml    = false;

    % Extract priors
    % Lambda ~ W(W0,nu0)
    % mu     ~ N(m0,inv(beta0*Lambda))
    nu0   = priors.nu;
    W0    = priors.W;
    beta0 = priors.beta;
    m0    = priors.m;
    if numel(nu0)   == 1, nu0   = repmat(  nu0,[1 K  ]); end
    if size(W0,3)   == 1, W0    = repmat(   W0,[1 1 K]); end
    if numel(beta0) == 1, beta0 = repmat(beta0,[1 K  ]); end
    if size(m0,2)   == 1, m0    = repmat(   m0,[1 K  ]); end
    nu    = zeros(size(nu0));
    W     = zeros(size(W0));
    beta  = zeros(size(beta0));
    m     = zeros(size(m0));
    out_L = zeros(4,K);

    if M >= min(nu0)+1,
        warning('SPM:Wishart','Can not normalise Wishart distribution (M=%d, nu_0=%f)', M,min(nu0));
    end
end

L      = -Inf;

for k=1:K,
    % Compute zeroeth moment
    s0 = 0;
    for i=1:numel(mom), s0 = s0 + mom(i).s0(k); end

    if ~ml,
        nu(k)    = nu0(k)   + s0;
        beta(k)  = beta0(k) + s0;
        W(:,:,k) = W0(:,:,k);
        m(:,k)   = m0(:,k);
    end

    L      = -Inf;
    for iter=1:1024,

        if ml,
            muk    = mu(:,k);
            Lambda = inv(Sigma(:,:,k));
        else
            muk    = m(:,k);
            Lambda = W(:,:,k)*nu(k);
        end

        s1     = zeros(M,1); % First moment
        S2     = zeros(M,M); % Second moment
        s1i    = zeros(M,1);
        S2i    = zeros(M,M);
        old_L  = L;

        % E[log q(h)]-------------
        L3     = 0;
        for i=1:numel(mom),
            if mom(i).s0(k),
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

        if ml,
            N(k)         = s0;
            mu(:,k)      = s1/s0;
            Sigma(:,:,k) = S2/s0 - mu(:,k)*mu(:,k)';
            L            = -0.5*s0*(logdet(Sigma(:,:,k)) + M*log(2*pi)) ...
                           -0.5*s0*trace((mu(:,k)*mu(:,k)' + S2/s0 - 2*mu(:,k)*s1'/s0)/Sigma(:,:,k)) ...
                           -L3;
           %fprintf('%d\t%g\n', iter,L);
        else
            m(:,k)   = (beta0(k)*m0(:,k) + s1)/beta(k);
            
            
            W(:,:,k) = inv(beta0(k)*m0(:,k)*m0(:,k)' + S2 - beta(k)*m(:,k)*m(:,k)' + inv(W0(:,:,k)));
            
            [V,D]    = eig(W(:,:,k));
            tol      = max(diag(D))*eps('single');
            D        = diag(max(diag(D),tol));
            W(:,:,k) = real(V*D*V');

%             
%             W(:,:,k) = inv(inv(W0(:,:,k))+ S2 ...
%                             - beta(k)\s1*s1'...
%                             + beta(k)\m0(:,k)*m0(:,k)'*beta0(k)*s0...
%                             - beta(k)\s1*m0(:,k)'*beta0(k)...
%                             - beta(k)\m0(:,k)*s1'*beta0(k));
%                         
            Eld = Elogdet(W(:,:,k),nu(k)); 
            if isempty(Eld)
                warning('W is not positive definite');
                fprintf(strcat('class #',num2str(k),'\n'));
                eig(W(:,:,k))
                beta(k)   = beta0(k);
                m(:,k)    = m0(:,k);
                nu(k)     = nu0(k);
                W0(:,:,k) = diag(eig(W0(:,:,k)))/2;
                W(:,:,k)  = W0(:,:,k);
                Eld       = Elogdet(W(:,:,k),nu(k));
                L1        =   0.5*s0*(Eld - M*log(2*pi) - trace((m(:,k)*m(:,k)' + inv(beta(k)*nu(k)*W(:,:,k)) + S2/s0 - 2*m(:,k)*s1'/s0)*nu(k)*W(:,:,k)));
                L2        =   0.5*M*(log(beta0(k)/(2*pi)) - beta0(k)/beta(k)) + 0.5*(nu0(k)-M)*Eld - logWishartDen(W0(:,:,k),nu0(k)) ...
                            - 0.5*nu(k)*beta0(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)) - 0.5*nu(k)*trace(W0(:,:,k)\W(:,:,k));
                L4        =   0.5*M*(log(beta(k)/(2*pi)) - 1 - nu(k)) + 0.5*(nu(k)-M)*Eld - logWishartDen(W(:,:,k),nu(k));
            end
            
%             % E[log p(g,h|mu,Lambda)]------------- 
            L1  = 0.5*s0*( Eld - M*log(2*pi) - ...
                trace( (m(:,k)*m(:,k)' + inv(beta(k)*nu(k)*W(:,:,k)) + S2/s0 - 2*m(:,k)*s1'/s0)*nu(k)*W(:,:,k) ) ); % (C.2)

%             L1 =   0.5*s0*(Eld - M*log(2*pi) - M/beta(k)) ...
%                  - 0.5*s0*(m(:,k)'*W(:,:,k)*m(:,k)) ...
%                  - 0.5*nu(k)*trace(W(:,:,k)*S2 - 2*s1*m(:,k)'*W(:,:,k)); % (C.2)
            
%             % E[log p(mu,Lambda)] -------------
            L2 = 0.5*M*(log(beta0(k)/(2*pi)) - beta0(k)/beta(k)) ...
                 + 0.5*(nu0(k)-M)*Eld - logWishartDen(W0(:,:,k),nu0(k)) ...
                 - 0.5*nu(k)*beta0(k)*(m(:,k)-m0(:,k))'*W(:,:,k)*(m(:,k)-m0(:,k)) ...
                 - 0.5*nu(k)*trace(W0(:,:,k)\W(:,:,k)); % (C.4)

%             L2 = 0.5*(M*log(beta0(k)/(2*pi)) - M*beta0(k)/beta(k)) ...
%                  + 2*logWishartDen(W0(:,:,k),nu0(k)) ...
%                  - 0.5*nu(k)*trace( (inv(W(:,:,k)) + beta0(k)*(m(:,k) - m0(:,k)*(m(:,k) - m0(:,k))'))*W(:,:,k) ) ...
%                  - Eld*(nu0(k) - M);
             
            % E[log q(mu,Lambda)] -------------
            L4 = 0.5*M*(log(beta(k)/(2*pi)) - 1 - nu(k)) + 0.5*(nu(k)-M)*Eld - logWishartDen(W(:,:,k),nu(k)); % (C.6) OK

            % Negative Free Energy
            L = L1 + L2 - L3 - L4;
           %fprintf('%d\t%g\t%g\t%g\t%g\t%g\n', iter,L1, L2, L3, L4, L);
        end
        if L-old_L < s0*M*1e-10, break; end
    end
    if ml,
        out_L(k)    = L;
    else
        out_L(:,k) = [L1 L2 L3 L4]';
    end
end


if ml,
    % Maximum likelihood estimates
    varargout{1} = N;
    varargout{2} = mu;
    varargout{3} = Sigma;
    varargout{4} = out_L;
else
    % Posterior distributions
    post.m       = m;
    post.W       = W;
    post.nu      = nu;
    post.beta    = beta;
    priors.nu    = nu0;
    priors.W     = W0;
    priors.beta  = beta0;
    priors.m     = m0;
    varargout{1} = post;
    varargout{2} = out_L;
    varargout{3} = priors;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = Elogdet(W,nu)
% E[log(det(Lambda))], where Lambda ~ W(W,nu)
M = size(W,1);
l = logdet(W); if isempty(l); e = []; return; end
e = l + M*log(2);
for m=1:M,
    e = e+psi((nu-m+1)/2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l,p] = logdet(W)
% log(det(W))
[cW,p] = chol(W); if p; l = []; return; end
l = 2*sum(log(diag(cW)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lZ = logWishartDen(W,nu)
% Log of normalising term of log W(W,nu)
M     = size(W,1);
if M >= nu+1,
   %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
    lZ = 0;
    return;
end
lGamM = M*(M-1)/4*log(pi);
for m=1:M,
    lGamM = lGamM+gammaln((nu+1-m)/2);
end
lZ    = 0.5*nu*(logdet(W) + M*log(2)) + lGamM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

