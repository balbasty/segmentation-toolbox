function lb = lowerbound(mom,po,pr)
mom = mom_John2Bishop(mom);  
s0  = mom.s0;
s1  = mom.s1;
S2  = mom.S2;

m = po.m;
n = po.n;
W = po.W;
b = po.b;

m0 = pr.m;
n0 = pr.n;
W0 = pr.W;
b0 = pr.b;

[N,K] = size(m0);

lb = 0;
for k=1:K
    logB0 = (n0(k)/2)*LogDet(inv(W0(:,:,k))) - (n0(k)*N/2)*log(2) ...
          - (N*(N-1)/4)*log(pi) - sum(gammaln(0.5*(n0(k)+1 -[1:N])));
  
    t1          = psi(0, 0.5*repmat(n(k)+1,N,1) - 0.5*[1:N]');
    logLamTilde = sum(t1) + N*log(2)  + LogDet(W(:,:,k));
    
    logBk = -(n(k)/2)*LogDet(W(:,:,k)) - (n(k)*N/2)*log(2)...
            - (N*(N-1)/4)*log(pi) - sum(gammaln(0.5*(n(k) + 1 - [1:N])));
    H     = -logBk - 0.5*(n(k) -N-1)*logLamTilde + 0.5*n(k)*N;

    trSW      = trace(n(k)*S2(:,:,k)*W(:,:,k));
    diff1     = s1(:,k) - m(:,k);
    xbarWxbar = diff1'*W(:,:,k)*diff1;

    diff1    = m(:,k) - m0(:,k);
    mWm      = b0(k)*n(k)*diff1'*W(:,:,k)*diff1; 
    trW0invW = trace(W0(:,:,k)\W(:,:,k));
    
    lb1 = 0.5*(s0(k).*(logLamTilde - N./b(k) - trSW - n(k).*xbarWxbar - N*log(2*pi)));
    lb2 = 0.5*(N*log(b0(k)/(2*pi)) + logLamTilde - N*(b0(k)./b(k)) - mWm);
    lb3 = logB0 + 0.5*((n0(k) - N - 1).*logLamTilde) - 0.5*(n(k).*trW0invW);    
    lb4 = 0.5*(logLamTilde + N.*log(b(k)/(2*pi))) - 0.5*N*K - H;
    
    lbk = lb1 + lb2 + lb3 - lb4;
    lb  = lb + lbk;
end
%=======================================================================