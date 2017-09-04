function [L,r,cp,Pf] = update_cp(f,bf,lnmu,cp,msk)

[Nf,D] = size(f);
K      = numel(cp.lnw);

f    = double(f);  % data
bf   = double(bf); % bias field
lnmu = double(lnmu); % template

msk1     = repmat(msk,1,D);
f(~msk1) = 0;

lnmu = reshape(lnmu,[Nf K]);

f = bf.*f;
f = f';

lnbias = log(prod(bf,2));
lnbias = repmat(lnbias,1,K);

% initialize variables
L = 0;

m0 = cp.pr.m;
b0 = cp.pr.b;
W0 = cp.pr.W;
n0 = cp.pr.n;

m = cp.po.m;
b = cp.po.b;
W = cp.po.W;
n = cp.po.n;

lnw = cp.lnw;

[lnPz,mgm] = multinomial(lnmu,lnw,msk);

% Calculate r
if nargout==4
    [r,~,~,Pf] = resps(f,K,Nf,D,b,n,m,W,lnbias,lnPz);
    return;
else
    [r,logr,logLambdaTilde] = resps(f,K,Nf,D,b,n,m,W,lnbias,lnPz);
end

ss = suffstats(f,r);

if nargout==3
    
    % Begin M step
    % compute new parameters
    [b,n,m,W,lnw] = vmstep(ss,K,b0,n0,m0,W0,mgm,lnw,cp.dow);

    cp.po.m = m;
    cp.po.b = b;
    cp.po.W = W;
    cp.po.n = n;
    cp.lnw  = lnw;
    
    lnPz = multinomial(lnmu,lnw,msk);
    
    [r,logr,logLambdaTilde] = resps(f,K,Nf,D,b,n,m,W,lnbias,lnPz);
    
    ss = suffstats(f,r);    
end

L = lowerbound(ss,K,D,r,logr,b0,n0,m0,W0,b,n,m,W,logLambdaTilde,lnbias,lnPz,msk);
%==========================================================================

%==========================================================================
function L = lowerbound(ss,K,D,r,logr,b0,n0,m0,W0,b,n,m,W,logLambdaTilde,lnbias,lnPz,msk)

Nk   = ss.s0;
xbar = ss.s1;
S    = ss.S2;

msk     = repmat(msk,1,K);
r(~msk) = 0;

% Various other parameters for different terms
H = 0;
logB0 = 0;
for k = 1:K
    W0inv = inv(W0(:,:,k));

    % B(lambda0)
    logB0 = logB0 + (n0(k)/2)*log(det(W0inv)) - (n0(k)*D/2)*log(2) ...
          - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n0(k)+1 -[1:D])));
  
  
    % sum(H(q(Lamba(k))))
    logBk = -(n(k)/2)*log(det(W(:,:,k))) - (n(k)*D/2)*log(2)...
            - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n(k) + 1 - [1:D])));
    H = H -logBk - 0.5*(n(k) -D-1)*logLambdaTilde(k) + 0.5*n(k)*D;
    % for Lt1 - third term
    trSW(k) = trace(n(k)*S(:,:,k)*W(:,:,k));
    diff = xbar(:,k) - m(:,k);
    xbarWxbar(k) = diff'*W(:,:,k)*diff;
    % for Lt4 - Fourth term
    diff = m(:,k) - m0(:,k);
    mWm(k) = b0(k)*n(k)*diff'*W(:,:,k)*diff; 
    trW0invW(k) = trace(W0inv*W(:,:,k));
end

Lt1 = 0.5*sum(Nk.*(logLambdaTilde - D./b...
    - trSW - n.*xbarWxbar - D*log(2*pi)));
Lt41 = 0.5*sum(D*log(b0/(2*pi)) + logLambdaTilde - D*sum(b0./b) - mWm);
Lt42 = logB0 + 0.5*sum((n0 - D - 1).*logLambdaTilde) - 0.5*sum(n.*trW0invW);
Lt4 = Lt41+Lt42;
Lt5 = sum(sum(r.*logr));
Lt7 = 0.5*sum(logLambdaTilde + D.*log(b/(2*pi))) - 0.5*D*K - H;
Lt3 = sum(sum(r.*lnbias));
Lt8 = sum(sum(r.*lnPz)); % (10.72)

%Bishop's Lower Bound
L = Lt1 + Lt3 + Lt4 - Lt5 - Lt7 + Lt8;
%==========================================================================
  
%==========================================================================
function [r,logr,logLambdaTilde,Pf] = resps(x,K,N,D,b,n,m,W,lnbias,lnPz)

E              = zeros(N,K);
logLambdaTilde = zeros(1,K);
for k = 1:K
    t1 = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
    logLambdaTilde(k) = sum(t1) + D*log(2)  + log(det(W(:,:,k)));
    
%     for nf = 1:N
%       % Calculate E
%       diff1   = x(:,nf) - m(:,k);
%       E(nf,k) = D/b(k) + n(k)*diff1'*W(:,:,k)*diff1;
%     end
    
    diff1  = bsxfun(@minus,x,m(:,k));
    Q      = chol(W(:,:,k))*diff1;
    E(:,k) = D/b(k) + n(k)*dot(Q,Q,1);
end
if nargout==4
    logRho     = repmat(0.5*logLambdaTilde, N,1) - 0.5*E  - D/2*log(2*pi);
    max_logRho = max(logRho,[],2);    
    Pf         = exp(bsxfun(@minus,logRho,max_logRho));
    logr       = 0;
    r          = 0;
else
    logRho    = repmat(0.5*logLambdaTilde, N,1) - 0.5*E  - D/2*log(2*pi) + lnbias + lnPz;
    logSumRho = logsumexp(logRho,2);
    logr      = logRho - repmat(logSumRho, 1,K);
    r         = exp(logr);
end
%==========================================================================

%==========================================================================
function [b,n,m,W,lnw] = vmstep(ss,K,b0,n0,m0,W0,mgm,lnw,dow)
D = size(m0,1);

Nk   = ss.s0;
xbar = ss.s1;
S    = ss.S2;

b = b0 + Nk;
n = n0 + Nk;
W = zeros(D,D,K);
for k = 1:K
    m(:,k) = (b0(k)*m0(:,k) + Nk(k).*xbar(:,k))./b(k);
    
    W0inv    = inv(W0(:,:,k));
    mlt1     = b0(k).*Nk(k)/(b0(k) + Nk(k));
    diff1    = xbar(:,k) - m0(:,k);
    W(:,:,k) = inv(W0inv + Nk(k)*S(:,:,k) + mlt1*(diff1*diff1'));
end  

if dow
    w   = (Nk + 1)./(mgm + K); % Bias the solution towards one
    w   = w/sum(w);
    lnw = log(w);
end
%==========================================================================

%==========================================================================
function [lnPz,mgm] = multinomial(lnmu,lnw,msk)
Pz = bsxfun(@plus,lnmu,lnw);

if nargout==2
    K = size(lnmu,2);        

    Pz1 = bsxfun(@times,exp(lnmu),exp(lnw));
    
    % Unified segmentation (27)
    mgm = 1./(sum(Pz1,2)); Pz1 = [];
    mgm = bsxfun(@times,lnmu,mgm);
    
    msk       = repmat(msk,1,K);
    mgm(~msk) = 0;
    mgm       = sum(mgm);
end

Pz   = softmax(Pz);
lnPz = log(Pz);
%==========================================================================

%==========================================================================
function mom = suffstats(X,Q,mom)
% Modified version of spm_suffstats: uses the same suffstats equations as
% in Bishop (10.51-53)
X = X';

if nargin<2 || isempty(Q),
    Q = ones(size(X,1),1);
end

K = size(Q,2);
M = size(X,2);
if M<=8,
    cast = @uint8;
    typ  = 'uint8';
elseif M<=16,
    cast = @uint16;
    typ  = 'uint16';
elseif M<=32,
    cast = @uint32;
    typ  = 'uint32';
elseif M<=64,
    cast = @uint64;
    typ  = 'uint64';
else,
    error('Too many dimensions.');
end

if nargin<3,
    % Create empty data structure
    mom = struct('ind',[],'s0',0,'s1',[],'S2',[]);
    for i=1:2^M,
        mom(i).ind = dec2bin(i-1,M)=='1'; % Indices
        Mi         = sum(mom(i).ind);
        mom(i).s0  = zeros(1,K);     % Zeroeth moments
        mom(i).s1  = zeros(Mi,K);    % First moments
        mom(i).S2  = zeros(Mi,Mi,K); % Second moments
    end
else
    % Check compatibility of data structure
    if (numel(mom)~=2^M) || (size(mom(1).s0,2)~=K),
        error('Incorrect moment dimensions');
    end
end
Q(isnan(Q))=0;
code = zeros([size(X,1),1],typ);
for i=1:M,
    % Mask for which combinations of data is missing
    code = bitor(code,bitshift(feval(cast,isfinite(X(:,i)) & (X(:,i)~=0)),(i-1)));
end
for i=2:numel(mom),
    msk0      = mom(i).ind;
    ind       = find(code==msk0*(2.^(0:(M-1))'));
    if ~isempty(ind),
        x  = X(ind,msk0);
        Nx = size(x,1);
        for k=1:K,
            q = Q(ind,k); 
            
            Nk             = sum(q);            
            mom(i).s0(1,k) = mom(i).s0(1,k) + Nk;
            
            xbar           = sum(repmat(q',M,1).*x',2)/Nk;                      
            mom(i).s1(:,k) = mom(i).s1(:,k) + xbar;                        
            
            diff1            = x' - repmat(xbar,1,Nx);
            diff2            = repmat(q',M,1).*diff1;    
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + (diff2*diff1')./Nk;    
        end
    end
end

mom = mom(end);
%==========================================================================

%==========================================================================
function s = logsumexp(b, dim)
% s = logsumexp(b) by Tom Minka
% Returns s(i) = log(sum(exp(b(:,i))))  while avoiding numerical underflow.
% s = logsumexp(b, dim) sums over dimension 'dim' instead of summing over rows

if nargin < 2 % if 2nd argument is missing
  dim = 1;
end

[B, junk] = max(b,[],dim);
dims = ones(1,ndims(b));
dims(dim) = size(b,dim);
b = b - repmat(B, dims);
s = B + log(sum(exp(b),dim));
i = find(~isfinite(B));
if ~isempty(i)
  s(i) = B(i);
end
%==========================================================================