function px = log_likelihoods(K,buf,mg,mog,cr)
if isempty(mg), mg = ones(1,K); end
if nargin<5,    cr = []; end

N = numel(buf.f);
I = numel(buf.code);
if isempty(cr)
    cr = zeros(I,N);
    for n1=1:N
        cr(buf.msk{n1},n1) = double(buf.f{n1}(:)).*double(buf.bf{n1}(:));
    end
end

if isfield(mog,'pr')
    vb = true;
    n  = mog.po.n;
    W  = mog.po.W;
    b  = mog.po.b;
    m  = mog.po.m;    
else
    vb = false;
    mn = mog.mn;
    vr = mog.vr;    
end

if vb
    nbf = ones([I N]);
    for n1=1:N
        nbf(buf.msk{n1},n1) = double(buf.bf{n1});
    end
    bf  = nbf; clear nbf
end

px = zeros(I,K);
for i=1:2^N    
    msk0 = dec2bin(i - 1,N)=='1';  
    ind  = find(buf.code==msk0*(2.^(0:(N-1))'));
    
    for k=1:K
        if vb                
            diff1     = bsxfun(@minus,cr(ind,msk0)',m(msk0,k));
            Q         = chol(W(msk0,msk0,k))*diff1;
            E         = N/b(k) + n(k)*dot(Q,Q,1);
            px(ind,k) = 0.5*(Elogdet(W(msk0,msk0,k),n(k)) - E') + log(mg(k)) - N/2*log(2*pi) + log(prod(bf(ind,msk0),2));
        else
            C         = chol(vr(msk0,msk0,k));
            diff1     = bsxfun(@minus,cr(ind,msk0),mn(msk0,k)')/C;
            px(ind,k) = log(mg(k)) - (N/2)*log(2*pi) - sum(log(diag(C))) - 0.5*sum(diff1.*diff1,2);
        end
    end
end
%=======================================================================