function L = log_likelihoods(f,bf,mg,gmm,code,cr)
if nargin<6, cr = []; end

K = numel(mg);
N = numel(f);
M = numel(f{1});

if isempty(cr)
    cr                 = zeros(M,N);
    for n=1:N, cr(:,n) = double(f{n}(:)).*double(bf{n}(:)); end
elseif iscell(cr)
    cr1                 = zeros(M,N);
    for n=1:N, cr1(:,n) = cr{n}; end    
    cr                  = cr1; clear cr1
end

if ~gmm.ml    
    nbf                 = ones([M N]);
    for n=1:N, nbf(:,n) = double(bf{n}); end
    bf                  = nbf; clear nbf
end

L = zeros(numel(f{1}),K);
% for n=2:2^N
%     msk0 = dec2bin(n - 1,N)=='1';
%     ind  = find(code==msk0*(2.^(0:(N - 1))'));
%     if ~isempty(ind)
        for k=1:K
            if gmm.ml
                C      = chol(gmm.vr(:,:,k));
                d      = bsxfun(@minus,cr,gmm.mn(:,k)')/C;
                L(:,k) = log(mg(k)) - (N/2)*log(2*pi) - sum(log(diag(C))) - 0.5*sum(d.*d,2);
            else
                d      = bsxfun(@minus,cr',gmm.po.m(:,k));
                Q      = chol(gmm.po.W(:,:,k))*d;
                E      = N/gmm.po.b(k) + gmm.po.n(k)*dot(Q,Q,1);
                L(:,k) = 0.5*(Elogdet(gmm.po.W(:,:,k),gmm.po.n(k)) - E') + log(mg(k)) - N/2*log(2*pi) + log(prod(bf,2));
            end
        end
    end
% end
%=======================================================================