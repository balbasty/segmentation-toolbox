function [lb,r,cp,lik] = update_cp(K,x,bf,lnmu,cp,msk,d) 

pr = cp.pr;
po = cp.po;

%--------------------------------------------------------------------------
% Mask

x    = double(x(msk));
bf   = double(bf(msk));
lnmu = double(lnmu);

[Nx,C] = size(x);

%--------------------------------------------------------------------------
% Apply bias field

x = bf(:).*x(:);

%--------------------------------------------------------------------------
% Initialize variables

lnLamTilde = zeros(1,K);
E          = zeros(Nx,K);
trSW       = zeros(1,K);
xbarWxbar  = zeros(1,K);
mWm        = zeros(1,K);
trW0invW   = zeros(1,K);

%--------------------------------------------------------------------------
% Set priors

m0    = pr.m;    
beta0 = pr.beta; 
W0    = pr.W;     
W0inv = zeros(C,C,K);
for k=1:K
    W0inv(:,:,k) = inv(W0(:,:,k));
end
nu0   = pr.nu;

%--------------------------------------------------------------------------
% Set posteriors

beta = po.beta;
nu   = po.nu;
m    = po.m;
W    = po.W;

%--------------------------------------------------------------------------
% Calculate responsibilities

[Pmu,mgm] = weighted_softmax(lnmu,log(cp.w),msk);

lnBias = log(prod(bf,2));
lnBias = repmat(lnBias,1,K);

for k=1:K
    t1            = psi(0,0.5*repmat(nu(k) + 1,C,1) - 0.5*(1:C)'); % (10.65)
    lnLamTilde(k) = sum(t1) + C*log(2)  + log(det(W(:,:,k)));      % (10.65)

    diff1  = bsxfun(@minus,x,m(:,k)');                                 % (10.64)
    E(:,k) = C/beta(k) + nu(k)*diff1.*(bsxfun(@times,W(:,:,k),diff1)); % (10.64)
end
lnRho = repmat(0.5*lnLamTilde,Nx,1) - 0.5*E - C/2*log(2*pi) + lnBias; % (10.46)

if nargout>3
    lb = 0;
    r  = 0;
    
    % VB-likelihood (10.67 w/o Pmu)
    lnSumRho = logsumexp(lnRho,2);
    lnr      = lnRho - repmat(lnSumRho, 1,K);
    msklik   = exp(lnr);
    
    lik = zeros([prod(d) K],'single');
    for k=1:K
        lik(msk,k) = single(msklik(:,k));
    end
    lik = reshape(lik,[d K]);
    
    return
end

% (10.46-49)
lnRho    = lnRho + log(Pmu); 
lnSumRho = logsumexp(lnRho,2);
lnr      = lnRho - repmat(lnSumRho, 1,K);
r        = exp(lnr);

if sum(~isfinite(r(:))),  warning('sum(~isfinite(r(:)))');  end

%--------------------------------------------------------------------------
% Compute sufficient statistics

s0 = sum(r); % (10.51)
s1 = zeros(C,K);
S2 = zeros(C,C,K);
for k=1:K
    s1(:,k)   = sum(bsxfun(@times,x,r(:,k)))/s0(k); % (10.52)
    diff1     = bsxfun(@minus,x,s1(:,k)');   % (10.53)
    S2(:,:,k) = sum(r(:,k).*diff1.^2)/s0(k); % (10.53)
end            

if nargout>2
    %----------------------------------------------------------------------
    % VM-step

    beta = beta0 + s0; % (10.60)
    nu   = nu0 + s0; % (10.63)
    m    = (repmat(beta0,C,1).*m0 + repmat(s0,C,1).*s1)./repmat(beta,C,1); % (10.61)
    for k=1:K
        mult1    = beta0(k).*s0(k)/(beta0(k) + s0(k));                         % (10.62)
        diff3    = s1(:,k) - m0(:,k);                                          % (10.62)
        W(:,:,k) = inv(W0inv(:,:,k) + s0(k)*S2(:,:,k) + mult1*(diff3*diff3')); % (10.62)
    end
    
    cp.po.m    = m;
    cp.po.W    = W;
    cp.po.nu   = nu;
    cp.po.beta = beta;

    %----------------------------------------------------------------------
    % Update weights     
    
    w = (s0 + 1)./(mgm + K); % Unified segmentation (27)
    w = w/sum(w);            % Unified segmentation (27)
    w = ones(size(w))/K;
    cp.w = w;
    
    Pmu = weighted_softmax(lnmu,log(cp.w),msk);
end

%--------------------------------------------------------------------------
% Compute lower bound

lnB0 = zeros(1,K);
H    = 0;
for k=1:K
    lnB0(k) = - (nu0(k)/2)*log(det(W0inv(:,:,k))) - (nu0(k)*C/2)*log(2) ...
              - (C*(C-1)/4)*log(pi) - sum(gammaln(0.5*(nu0(k) + 1 - (1:C)))); % (B.79)
  
    lnBk = - (nu(k)/2)*log(det(W(:,:,k))) - (nu(k)*C/2)*log(2)...
           - (C*(C-1)/4)*log(pi) - sum(gammaln(0.5*(nu(k) + 1 - (1:C)))); % (B.79)
    H    = H - lnBk - 0.5*(nu(k) - C - 1)*lnLamTilde(k) + 0.5*nu(k)*C;       % (B.82)
    
    trSW(k) = trace(nu(k)*S2(:,:,k)*W(:,:,k)); % (10.71)
   
    diff1        = s1(:,k) - m(:,k);      % (10.71)
    xbarWxbar(k) = diff1'*W(:,:,k)*diff1; % (10.71)
   
    diff1       = m(:,k) - m0(:,k);             % (10.74)
    mWm(k)      = diff1'*W(:,:,k)*diff1;        % (10.74)
    trW0invW(k) = trace(W0inv(:,:,k)*W(:,:,k)); % (10.74)
end

lb1  = 0.5*sum(s0.*(lnLamTilde - C./beta - trSW - nu.*xbarWxbar - C*log(2*pi))); % (10.71)
lb2  = sum(sum(r.*log(Pmu),2),1); % (10.72)
lb31 = 0.5*sum(C*log(beta0/(2*pi)) + lnLamTilde - C*beta0./beta - beta0.*nu.*mWm); % (10.74)
lb32 = sum(lnB0) + 0.5*sum((nu0 - C - 1).*lnLamTilde) - 0.5*sum(nu.*trW0invW);     % (10.74)
lb3  = lb31 + lb32;                                                                % (10.74)
lb4  = sum(sum(r.*lnBias,2),1);
lb5  = sum(sum(r.*lnr,2),1); % (10.75)
lb6  = 0.5*sum(lnLamTilde + C.*log(beta/(2*pi))) - 0.5*C*K - H; % (10.77)

if ~isfinite(lb1),  warning('~isfinite(lb1)');  end
if ~isfinite(lb2),  warning('~isfinite(lb2)');  end
if ~isfinite(lb31), warning('~isfinite(lb31)'); end
if ~isfinite(lb32), warning('~isfinite(lb32)'); end
if ~isfinite(lb4),  warning('~isfinite(lb4)');  end
if ~isfinite(lb5),  warning('~isfinite(lb5)');  end
if ~isfinite(lb6),  warning('~isfinite(lb6)');  end

lb = lb1 + lb2 + lb3 + lb4 - lb5 - lb6;

%--------------------------------------------------------------------------

or = r;
r  = ones([prod(d) K],'single')/K;
for k=1:K
    r(msk,k) = single(or(:,k));
end
%==========================================================================

%==========================================================================
function s = logsumexp(b, dim)
if nargin < 2
  dim = 1;
end

[B,~] = max(b,[],dim);
dims = ones(1,ndims(b));
dims(dim) = size(b,dim);
b = b - repmat(B, dims);
s = B + log(sum(exp(b),dim));
i = find(~isfinite(B));
if ~isempty(i)
  s(i) = B(i);
end
%==========================================================================