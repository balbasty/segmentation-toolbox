function [L,r,cp,Pf] = update_cp(f,bf,lnmu,cp)

missingdata = true;

K = numel(cp.lnw);

lnmu = reshape(lnmu,[size(f,1) K]);

f = bf.*f;

lnbf = log(prod(bf,2));
lnbf = repmat(lnbf,1,K);

lnPz = multinomial(lnmu,cp.lnw);

% Calculate r
if nargout==4
    if missingdata
        [r,~,Pf] = resps_missing(f,cp.po,lnbf,lnPz);
    else
        [r,~,Pf] = resps(f,cp.po,lnbf,lnPz);
    end
    L = 0;
    return;
else
    if missingdata
        [r,lnr] = resps_missing(f,cp.po,lnbf,lnPz);
    else
        [r,lnr] = resps(f,cp.po,lnbf,lnPz);    
    end
end

mom = suffstats(f,r);
if ~missingdata
    mom = mom(end);
end

if missingdata
    msk = sum(f>0,2)>0;
else
    msk = sum(f>0,2)==size(f,2);
end

if nargout==3
    
    % Begin M step
    % compute new parameters
    if missingdata
        [~,cp.po] = VBGaussiansFromSuffStats(mom,cp.pr,cp.po);
    else
        cp.po     = vmstep(mom,cp.pr);   
    end        
    
    if cp.dow       
        b = zeros([nnz(msk) K],'single');
        for k=1:K
            b(:,k) = exp(lnmu(msk,k));
        end

        w  = exp(cp.lnw);
        bw = bsxfun(@times,b,w);    
        bw = 1./(sum(bw,2));
        
        mgm = bsxfun(@times,b,bw); bw = [];
        mgm = sum(mgm);
        
        s0 = 0;
        for m=1:numel(mom)
            s0 = s0 + mom(m).s0;
        end
    
        w      = (s0 + 1)./(mgm + K);
        w      = w/sum(w);        
        cp.lnw = log(w);
        
        if sum(~isfinite(w))
           error('sum(~isfinite(w))');
        end
        
        lnPz = multinomial(lnmu,cp.lnw);
    end
    
    if missingdata
        [r,lnr] = resps_missing(f,cp.po,lnbf,lnPz);
    else
        [r,lnr] = resps(f,cp.po,lnbf,lnPz);
    end
    
    mom = suffstats(f,r); 
    if ~missingdata
        mom = mom(end);
    end
end

    if missingdata
        L = VBGaussiansFromSuffStats(mom,cp.pr,cp.po);
    else
        L = lowerbound(mom,cp.pr,cp.po);
    end
    
if ~missingdata
    msk     = repmat(msk,1,K);
    r(~msk) = NaN;
end

L5 = nansum(nansum(r.*lnr));
L6 = nansum(nansum(r.*lnbf));
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
function [L,po] = VBGaussiansFromSuffStats(mom,pr,po,verbose)
if nargin<4, verbose=false; end

K = size(mom(1).s0,2);
D = numel(mom(1).ind);

% Get priors
m0 = pr.m;
b0 = pr.b;
W0 = pr.W;
n0 = pr.n;

% Get posteriors
n = po.n;
W = po.W;
b = po.b;
m = po.m;

if verbose
    fprintf('----------------------\n');
end

L  = 0;
for k=1:K,
    
    s0 = 0;
    for i=1:numel(mom), 
        s0 = s0 + mom(i).s0(k); 
    end           
    
    Lk = -Inf;
    for iter=1:1024,
        oLk = Lk; 
        
        mu = m(:,k);
        P  = W(:,:,k)*n(k);
            
        s1  = zeros([D,1],'single');
        S2  = zeros(D,'single');        
        s1i = zeros([D,1],'single');
        S2i = zeros([D,D],'single');         
        L5  = 0;
        for i=1:numel(mom),
            
            if mom(i).s0(k),
                
                ind = mom(i).ind;
                s0m = mom(i).s0(k);
                s1m = mom(i).s1(:,k);
                S2m = mom(i).S2(:,:,k);
                        
                [s1m,S2m] = mom_Bishop2John(s0m,s1m,S2m);
                
                mux            = mu( ind,:);
                muy            = mu(~ind,:);
                Pyy            = P(~ind,~ind);
                Pyx            = P(~ind, ind);
                R              = Pyy\Pyx;
                Ex             = s1m/s0m;
                Exx            = S2m/s0m;
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = s1m;
                S2i( ind, ind) = S2m;
                s1i(~ind)      = s0m*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = s0m*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Pyy));
                S2i( ind,~ind) = s0m*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;
                
                L5 = L5 + 0.5*s0m*(logdet(Pyy) - size(Pyy,1)*(1+log(2*pi)));
            end
        end

        [s1,S2] = mom_John2Bishop(s0,s1,S2);
        
        W0inv = inv(W0(:,:,k));
        
        if nargout>1
            % VM-step----------------------------------------------------------                
            b(k) = b0(k) + s0;

            n(k) = n0(k) + s0;

            m(:,k) = (b0(k)*m0(:,k) + s0.*s1)./b(k);
            
            mlt1     = b0(k).*s0/(b0(k) + s0);
            diff1    = s1 - m0(:,k);
            W(:,:,k) = inv(W0inv + s0*S2 + mlt1*(diff1*diff1'));
            
%             [V,D1]   = eig(W(:,:,k));
%             tol      = max(diag(D1))*eps('single');
%             D1       = diag(max(diag(D1),tol));
%             W(:,:,k) = real(V*D1*V');
        end
        
        % Compute objective function---------------------------------------        
        logB0 = (n0(k)/2)*logdet(W0inv) - (n0(k)*D/2)*log(2) ...
              - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n0(k)+1 -[1:D])));

        t1          = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
        logLamTilde = sum(t1) + D*log(2)  + logdet(W(:,:,k));

        logBk = -(n(k)/2)*logdet(W(:,:,k)) - (n(k)*D/2)*log(2)...
                - (D*(D-1)/4)*log(pi) - sum(gammaln(0.5*(n(k) + 1 - [1:D])));
        H     = -logBk - 0.5*(n(k) -D-1)*logLamTilde + 0.5*n(k)*D;

        trSW      = trace(n(k)*S2*W(:,:,k));
        diff      = s1 - m(:,k);
        xbarWxbar = diff'*W(:,:,k)*diff;

        diff     = m(:,k) - m0(:,k);
        mWm      = b0(k)*n(k)*diff'*W(:,:,k)*diff; 
        trW0invW = trace(W0inv*W(:,:,k));

        L1 = 0.5*(s0.*(logLamTilde - D./b(k) - trSW - n(k).*xbarWxbar - D*log(2*pi)));
        L2 = 0.5*(D*log(b0(k)/(2*pi)) + logLamTilde - D*(b0(k)./b(k)) - mWm);
        L3 = logB0 + 0.5*((n0(k) - D - 1).*logLamTilde) - 0.5*(n(k).*trW0invW);    
        L4 = 0.5*(logLamTilde + D.*log(b(k)/(2*pi))) - 0.5*D*K - H;

        Lk = L1 + L2 + L3 - L4 - L5;        

        if verbose
            fprintf('%d\t%g\n', iter,Lk);
        end
        if abs(Lk-oLk) < 1e-12, 
            L = L + Lk;
            break; 
        end
    end
end

% update posteriors
po.m = m;
po.n = n;
po.W = W;
po.b = b;
%==========================================================================

%==========================================================================
function [r,lnr,Pf] = resps(f,po,lnbf,lnPz)
m = po.m;
n = po.n;
W = po.W;
b = po.b;

[D,K] = size(m);
N     = size(lnPz,1);

E          = zeros([N,K],'single');
lnLamTilde = zeros([1,K],'single');
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
    lnRho     = repmat(0.5*lnLamTilde, N,1) - 0.5*E  - D/2*log(2*pi) + lnbf + lnPz;
    lnSumRho  = logsumexp(lnRho,2);
    lnr       = lnRho - repmat(lnSumRho, 1,K);
    r         = exp(lnr);
end
%==========================================================================

%==========================================================================
function [r,lnr,Pf] = resps_missing(f,po,lnbf,lnPz)
m = po.m;
n = po.n;
W = po.W;
b = po.b;

[D,K] = size(m);

if D<=8,
    cast = @uint8;
    typ  = 'uint8';
elseif D<=16,
    cast = @uint16;
    typ  = 'uint16';
elseif D<=32,
    cast = @uint32;
    typ  = 'uint32';
elseif D<=64,
    cast = @uint64;
    typ  = 'uint64';
else,
    error('Too many dimensions.');
end

code = zeros([size(f,1),1],typ);
for i=1:D,
    code = bitor(code,bitshift(feval(cast,isfinite(f(:,i)) & (f(:,i)~=0)),(i-1)));
end

lnRho = NaN([size(f,1),K],'single');
for i=1:2^D % For combinations of missing data
    
    msk0 = dec2bin(i-1,D)=='1';
    ind  = find(code==msk0*(2.^(0:(D-1))'));      
    
    if ~isempty(ind)
        fi  = f(ind,msk0);       
        Nfi = size(fi,1);
        
        E          = zeros([Nfi,K],'single');
        lnLamTilde = zeros([1,K],'single');
        for k=1:K % For classes
           
            t1            = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
            lnLamTilde(k) = sum(t1) + D*log(2) + logdet(W(msk0,msk0,k));

            diff1  = bsxfun(@minus,fi',m(msk0,k));
            Q      = chol(W(msk0,msk0,k))*diff1;
            E(:,k) = D/b(k) + n(k)*dot(Q,Q,1)';
        end
        
        if nargout==3
            lnRho(ind,:) = repmat(0.5*lnLamTilde,Nfi,1) - 0.5*E  - D/2*log(2*pi);
        else
            lnRho(ind,:) = repmat(0.5*lnLamTilde,Nfi,1) - 0.5*E  - D/2*log(2*pi) + lnPz(ind,:) + lnbf(ind,:);
        end
    end
end

if nargout==3
    max_lnRho = nanmax(lnRho,[],2);    
    Pf        = exp(bsxfun(@minus,lnRho,max_lnRho));
    lnr       = 0;
    r         = 0;    
%     lnSumRho  = logsumexp(lnRho,2);
%     lnr       = lnRho - repmat(lnSumRho,1,K);
%     Pf        = exp(lnr);    
%     r         = 0;    
else
    lnSumRho  = logsumexp(lnRho,2);
    lnr       = lnRho - repmat(lnSumRho,1,K);
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
W = zeros([D,D,K],'single');
for k = 1:K
    m(:,k) = (b0(k)*m0(:,k) + s0(k).*s1(:,k))./b(k);
    
    W0inv    = inv(W0(:,:,k));
    mlt1     = b0(k).*s0(k)/(b0(k) + s0(k));
    diff1    = s1(:,k) - m0(:,k);
    W(:,:,k) = inv(W0inv + s0(k)*S2(:,:,k) + mlt1*(diff1*diff1'));
    
    [V,D]    = eig(W(:,:,k));
    tol      = max(diag(D))*eps('single');
    D        = diag(max(diag(D),tol));
    W(:,:,k) = real(V*D*V');
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
        mom(i).s0  = zeros([1,K],'single');     % Zeroeth moments
        mom(i).s1  = zeros([Mi,K],'single');    % First moments
        mom(i).S2  = zeros([Mi,Mi,K],'single'); % Second moments
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
        x = X(ind,msk0);
        for k=1:K,
            q = Q(ind,k); 
            
            s0             = sum(q) + eps;            
            mom(i).s0(1,k) = mom(i).s0(1,k) + s0;
            
            if s0==0
                error('s0==0');
            end
            
            s1             = (sum(bsxfun(@times,q,x))/s0)';                      
            mom(i).s1(:,k) = mom(i).s1(:,k) + s1;                        
            
            diff1            = bsxfun(@minus,x,s1');
            diff2            = bsxfun(@times,q,diff1);
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + (diff2'*diff1)./s0;    
        end
    end
end
%==========================================================================

%==========================================================================
function s = logsumexp(b, dim)
B         = nanmax(b,[],dim);
dims      = ones(1,ndims(b));
dims(dim) = size(b,dim);
b         = b - repmat(B, dims);
s         = B + log(nansum(exp(b),dim));
%==========================================================================

%==========================================================================
function [s1B,S2B] = mom_John2Bishop(s0,s1J,S2J)
s1B = s1J/s0;
S2B = S2J/s0 - (s1J/s0)*(s1J/s0)';
%==========================================================================

%==========================================================================
function [s1J,S2J] = mom_Bishop2John(s0,s1B,S2B)
s1J = s0*s1B;
S2J = s0*S2B + s0*(s1J/s0)*(s1J/s0)';
%==========================================================================