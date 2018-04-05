function L = log_likelihoods(f,bf,mg,gmm,msk,code,lkp,cr)
if nargin<8, cr = []; end

K = numel(mg);
N = numel(f);
M = numel(code);

% Compute Bx
%--------------------------------------------------------------------------
if isempty(cr)
    cr                       = NaN(M,N);
    for n=1:N, cr(msk{n},n)  = double(f{n}(:)).*double(bf{n}(:)); end
elseif iscell(cr)
    cr1                      = NaN(M,N);
    for n=1:N, cr1(msk{n},n) = cr{n}; end    
    cr                       = cr1; clear cr1
end

% Compute log|B|
%--------------------------------------------------------------------------
nbf                      = NaN([M N]);
for n=1:N, nbf(msk{n},n) = double(bf{n}); end

% Compute likelihoods
%--------------------------------------------------------------------------
L = NaN(M,K);
for n=2:2^N
    msk0                 = dec2bin(n - 1,N)=='1';
    ind                  = find(code==msk0*(2.^(0:(N - 1))'));
    if ~isempty(ind)
        for k=1:K
            if any(lkp.rem==lkp.part(k))
                L(ind,k) = -Inf;
            else
                d        = bsxfun(@minus,cr(ind,msk0)',gmm.po.m(msk0,k));
                Q        = chol(gmm.po.W(msk0,msk0,k))*d;
                E        = N/gmm.po.b(k) + gmm.po.n(k)*dot(Q,Q,1);
                L(ind,k) = 0.5*(Elogdet(gmm.po.W(msk0,msk0,k),gmm.po.n(k)) - E') + log(mg(k)) - N/2*log(2*pi) + log(prod(nbf(ind,msk0),2));
            end
        end
    end
end
%=======================================================================