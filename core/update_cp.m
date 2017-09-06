function [L,r,cp,Pf] = update_cp(f,bf,lnmu,cp,msk)

[Nf,D] = size(f);
K      = numel(cp.lnw);

f    = double(f);  % data
bf   = double(bf); % bias field
lnmu = double(lnmu); % template
lnmu = reshape(lnmu,[Nf K]);

% OBS: msk!
% msk     = sum(f,2)>0;
% f(~msk) = 0;

f = bf.*f;

lnbias = log(prod(bf,2));
lnbias = repmat(lnbias,1,K);

lnPz = multinomial(lnmu,cp.lnw);

% Calculate r
if nargout==4
    [r,~,Pf] = resps(f,cp.po,lnbias,lnPz);
    L        = 0;
    return;
else
    [r,logr] = resps(f,cp.po,lnbias,lnPz);
end

mom = suffstats(f,r);

if nargout==3
    
    if cp.dow       
        b = zeros([nnz(msk) K]);
        for k=1:K
            b(:,k) = exp(lnmu(msk,k));
        end

        bw = bsxfun(@times,b,exp(cp.lnw));    
        bw = 1./(sum(bw,2));
        
        mgm = bsxfun(@times,b,bw); bw = [];
        mgm = sum(mgm);
        
        w      = (mom.s0 + 1)./(mgm + K); % Bias the solution towards one
        w      = w/sum(w);
        cp.lnw = log(w);
        
        lnPz = multinomial(lnmu,cp.lnw);
    end

    % Begin M step
    % compute new parameters
    cp.po = vmstep(mom,cp.pr);   
    
    [r,logr] = resps(f,cp.po,lnbias,lnPz);
    
    mom = suffstats(f,r);    
end

L = lowerbound(mom,cp.pr,cp.po);

msk     = repmat(msk,1,K);
r(~msk) = NaN;

L5 = nansum(nansum(r.*logr));
L6 = nansum(nansum(r.*lnbias));
L7 = nansum(nansum(r.*lnPz));

%Bishop's Lower Bound
L = L - L5 + L6 + L7;
%==========================================================================

%==========================================================================
function L = lowerbound(mom,pr,po)

s0 = mom.s0;
s1 = mom.s1;
S2 = mom.S2;

m = po.m;
n = po.n;
W = po.W;
b = po.b;

m0 = pr.m;
n0 = pr.n;
W0 = pr.W;
b0 = pr.b;

[D,K] = size(m0);

L = 0;
for k = 1:K
    W0inv = inv(W0(:,:,k));

    logB0 = (n0(k)/2)*logdet(W0inv) - (n0(k)*D/2)*log(2) ...
          - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n0(k)+1 -[1:D])));
  
    t1          = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
    logLamTilde = sum(t1) + D*log(2)  + logdet(W(:,:,k));
    
    logBk = -(n(k)/2)*logdet(W(:,:,k)) - (n(k)*D/2)*log(2)...
            - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n(k) + 1 - [1:D])));
    H     = -logBk - 0.5*(n(k) -D-1)*logLamTilde + 0.5*n(k)*D;

    trSW      = trace(n(k)*S2(:,:,k)*W(:,:,k));
    diff      = s1(:,k) - m(:,k);
    xbarWxbar = diff'*W(:,:,k)*diff;

    diff     = m(:,k) - m0(:,k);
    mWm      = b0(k)*n(k)*diff'*W(:,:,k)*diff; 
    trW0invW = trace(W0inv*W(:,:,k));
    
    L1 = 0.5*(s0(k).*(logLamTilde - D./b(k) - trSW - n(k).*xbarWxbar - D*log(2*pi)));
    L2 = 0.5*(D*log(b0(k)/(2*pi)) + logLamTilde - D*(b0(k)./b(k)) - mWm);
    L3 = logB0 + 0.5*((n0(k) - D - 1).*logLamTilde) - 0.5*(n(k).*trW0invW);    
    L4 = 0.5*(logLamTilde + D.*log(b(k)/(2*pi))) - 0.5*D*K - H;
    
    Lk = L1 + L2 + L3 - L4;
    L  = L + Lk;
end
%==========================================================================
  
%==========================================================================
function [r,lnr,Pf] = resps(f,po,lnbias,lnPz)

m = po.m;
n = po.n;
W = po.W;
b = po.b;

[D,K] = size(m);
N     = size(lnPz,1);

E              = zeros(N,K);
lnLamTilde = zeros(1,K);
for k = 1:K
    t1 = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
    lnLamTilde(k) = sum(t1) + D*log(2)  + logdet(W(:,:,k));    
    
    diff1  = bsxfun(@minus,f',m(:,k));
    Q      = chol(W(:,:,k))*diff1;
    E(:,k) = D/b(k) + n(k)*dot(Q,Q,1);
end
if nargout==3
    lnRho     = repmat(0.5*lnLamTilde, N,1) - 0.5*E  - D/2*log(2*pi);
    max_lnRho = max(lnRho,[],2);    
    Pf        = exp(bsxfun(@minus,lnRho,max_lnRho));
    lnr       = 0;
    r         = 0;
else
    lnRho     = repmat(0.5*lnLamTilde, N,1) - 0.5*E  - D/2*log(2*pi) + lnbias + lnPz;
    lnSumRho  = logsumexp(lnRho,2);
    lnr       = lnRho - repmat(lnSumRho, 1,K);
    r         = exp(lnr);
end
%==========================================================================

%==========================================================================
function po = vmstep(mom,pr)

s0 = mom.s0;
s1 = mom.s1;
S2 = mom.S2;

m0 = pr.m;
n0 = pr.n;
W0 = pr.W;
b0 = pr.b;

[D,K] = size(m0);

b = b0 + s0;
n = n0 + s0;
W = zeros(D,D,K);
for k = 1:K
    m(:,k) = (b0(k)*m0(:,k) + s0(k).*s1(:,k))./b(k);
    
    W0inv    = inv(W0(:,:,k));
    mlt1     = b0(k).*s0(k)/(b0(k) + s0(k));
    diff1    = s1(:,k) - m0(:,k);
    W(:,:,k) = inv(W0inv + s0(k)*S2(:,:,k) + mlt1*(diff1*diff1'));
end  

po.m = m;
po.b = b;
po.W = W;
po.n = n;
%==========================================================================

%==========================================================================
function lnPz = multinomial(lnmu,lnw)
Pz   = bsxfun(@plus,lnmu,lnw);
Pz   = softmax(Pz);
lnPz = log(Pz);
%==========================================================================

%==========================================================================
function mom = suffstats(X,Q,mom)
% Modified version of spm_suffstats: uses the same suffstats equations as
% in Bishop (10.51-53)
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
            
            xbar           = (sum(bsxfun(@times,q,x))/Nk)';                      
            mom(i).s1(:,k) = mom(i).s1(:,k) + xbar;                        
            
            diff1            = bsxfun(@minus,x,xbar');
            diff2            = bsxfun(@times,q,diff1);
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + (diff2'*diff1)./Nk;    
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