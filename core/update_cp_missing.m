function [L,r,cp,Pf] = update_cp_missing(f,bf,lnmu,cp,msk)
K = numel(cp.lnw);

f    = double(f);    % data
bf   = double(bf);   % bias field
lnmu = double(lnmu); % log of template
lnmu = reshape(lnmu,[size(f,1) K]);

% msk1     = repmat(msk,1,D);
% f(~msk1) = 0;

f = bf.*f;

lnbf = log(prod(bf,2));
lnbf = repmat(lnbf,1,K);

% initialize variables
lnw = cp.lnw;
pr  = cp.pr;
po  = cp.po;
dow = cp.dow;

lnPz = multinomial(lnmu,lnw);

if nargout==4
    [r,~,Pf] = resps(f,K,po,lnbf,lnPz);
    L        = 0;
    return;
else
    [r,lnr] = resps(f,K,po,lnbf,lnPz);
end

mom = suffstats(f,r);

if nargout == 3
    [~,po,mg] = VBGaussiansFromSuffStats(mom,pr,po);

    if dow
        msk = isfinite(r(:,1));
%         msk = sum(f,2)>0;

        b = zeros([nnz(msk) K]);
        for k=1:K
            b(:,k) = exp(lnmu(msk,k));
        end
        
        w  = exp(lnw);
        bw = bsxfun(@times,b,w);    
        bw = 1./(sum(bw,2));
        
        mgm = bsxfun(@times,b,bw); bw = [];
        mgm = sum(mgm,1);

        w   = (mg + mgm)./(mgm); % Bias the solution towards one
%         w   = (mg + 1)./(mgm + K); 
%         w   = w/sum(w);
        lnw = log(w);
        
        lnPz = multinomial(lnmu,lnw);

        [r,lnr] = resps(f,K,po,lnbf,lnPz);

        mom = suffstats(f,r);
    end

    cp.po.m = po.m;
    cp.po.b = po.b;
    cp.po.W = po.W;
    cp.po.n = po.n;
    cp.lnw  = lnw;
end

L = VBGaussiansFromSuffStats(mom,pr,po);

L = L - nansum(nansum(r.*lnr));
L = L + nansum(nansum(r.*lnbf));
L = L + nansum(nansum(r.*lnPz));
%==========================================================================
  
%==========================================================================
function [L,po,mg] = VBGaussiansFromSuffStats(mom,pr,po)
K = size(mom(1).s0,2);
D = numel(mom(1).ind);

% Get priors
m0 = pr.m;
b0 = pr.b;
W0 = pr.W;
n0 = pr.n;

% Initialise posteriors
n = po.n;
W = po.W;
b = po.b;
m = po.m;

mg = zeros(1,K);
L  = 0;
for k=1:K,
    
    s0 = 0;
    for i=1:numel(mom), 
        s0 = s0 + mom(i).s0(k); 
    end            
    mg(k) = s0;
    
    Lk = -Inf;
    for iter=1:1024,
        oLk = Lk; 
        
%         if iter>1
            mu = m(:,k);
            P  = (W(:,:,k))*n(k);
%         else
%             P  = inv(eye(D));
%             mu = zeros(D,1);     
%         end
                
        s1  = zeros(D,1);
        S2  = zeros(D);        
        s1i = zeros(D,1);
        S2i = zeros(D,D);         
        L5  = 0;
        for i=1:numel(mom),
            
            if mom(i).s0(k),
                ind            = mom(i).ind;
                mux            = mu( ind,:);
                muy            = mu(~ind,:);
                Pyy            = P(~ind,~ind);
                Pyx            = P(~ind, ind);
                R              = Pyy\Pyx;
                Ex             = mom(i).s1(:,k)/mom(i).s0(k);
                Exx            = mom(i).S2(:,:,k)/mom(i).s0(k);
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = mom(i).s1(:,k);
                S2i( ind, ind) = mom(i).S2(:,:,k);
                s1i(~ind)      = mom(i).s0(k)*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = mom(i).s0(k)*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Pyy));
                S2i( ind,~ind) = mom(i).s0(k)*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;
                
                L5 = L5 + 0.5*mom(i).s0(k)*(logdet(Pyy) - size(Pyy,1)*(1+log(2*pi)));
            end
        end

        W0inv = inv(W0(:,:,k));
        
        if nargout>1
            % VM-step----------------------------------------------------------                
            b(k) = b0(k) + s0;

            n(k) = n0(k) + s0;

            m(:,k) = (b0(k)*m0(:,k) + s0.*s1)./b(k);
            
            mlt1     = b0(k).*s0/(b0(k) + s0);
            diff1    = s1 - m0(:,k);
            W(:,:,k) = inv(W0inv + s0*S2 + mlt1*(diff1*diff1'));
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

%         fprintf('%d\t%g\n', iter,Lk);
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
function [R,lnr,Pf] = resps(F,K,po,lnbf,lnPz)
D = size(F,2);

m = po.m;
n = po.n;
W = po.W;
b = po.b;

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

code = zeros([size(F,1),1],typ);
for i=1:D,
    code = bitor(code,bitshift(feval(cast,isfinite(F(:,i)) & (F(:,i)~=0)),(i-1)));
end

logRho = zeros(size(F,1),K);
for i=2:2^D % For combinations of missing data
    
    msk0 = dec2bin(i-1,D)=='1';
    ind  = find(code==msk0*(2.^(0:(D-1))'));  
    
    if ~msk0
        msk0 = dec2bin(2^D-1,D)=='1'; % Is this an okay approach?
    end
    
    if ~isempty(ind)
        Fi = F(ind,msk0);       
        Ni = size(Fi,1);
        
        E           = zeros(Ni,K);
        logLamTilde = zeros(1,K);
        for k=1:K % For classes
           
            t1             = psi(0, 0.5*repmat(n(k)+1,D,1) - 0.5*[1:D]');
            logLamTilde(k) = sum(t1) + D*log(2) + logdet(W(msk0,msk0,k));

            diff1  = bsxfun(@minus,Fi',m(msk0,k));
            Q      = chol(W(msk0,msk0,k))*diff1;
            E(:,k) = D/b(k) + n(k)*dot(Q,Q,1)';
        end
        
        if nargout==4
            logRho(ind,:) = repmat(0.5*logLamTilde,Ni,1) - 0.5*E  - D/2*log(2*pi);
        else
            logRho(ind,:) = repmat(0.5*logLamTilde,Ni,1) - 0.5*E  - D/2*log(2*pi) + lnPz(ind,:) + lnbf(ind,:);
        end
    end
end

if nargout==3
    max_logRho        = nanmax(logRho,[],2);    
    Pf                = exp(bsxfun(@minus,logRho,max_logRho));
    lnr              = 0;
    R                 = 0;
else
    logSumRho         = logsumexp(logRho,2);
    lnr              = logRho - repmat(logSumRho,1,K);
    R                 = exp(lnr);
end
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
%==========================================================================

%==========================================================================
function s = logsumexp(b, dim)
% s = logsumexp(b) by Tom Minka
% Returns s(i) = log(sum(exp(b(:,i))))  while avoiding numerical underflow.
% s = logsumexp(b, dim) sums over dimension 'dim' instead of summing over rows

if nargin < 2 % if 2nd argument is missing
  dim = 1;
end

B = nanmax(b,[],dim);
dims = ones(1,ndims(b));
dims(dim) = size(b,dim);
b = b - repmat(B, dims);
s = B + log(nansum(exp(b),dim));
i = find(~isfinite(B));
if ~isempty(i)
  s(i) = B(i);
end
%==========================================================================