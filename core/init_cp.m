function cp = init_cp(pthf,pthmu,mu,K,d,do)
[N,D] = size(pthf);

wpo   = zeros(1,K,N);
mupo  = zeros(D,K,N);
Sigpo = zeros(D,D,K,N);
Nf    = zeros(N,1);

for n1=1:N
    f   = zeros([prod(d) D]);
    msk = zeros([prod(d),D],'logical');
    for c=1:D
        Nii    = nifti(pthf{n1,c});
        f(:,c) = reshape(Nii.dat(:,:,:),[],1);
        
        msk(:,c) = get_msk(f(:,c));
    end
    Nf(n1) = nnz(f(:))/D;    

    msk     = sum(msk,2)==D;
    f(~msk) = 0;
    
    [mu,Sig,w] = Kmeans(f,K);

    % To ensure positive semi-definite covariance matrix
    for k=1:K
        B          = (Sig(:,:,k) + Sig(:,:,k)')/2; 
        [U,Sigma]  = eig(B); 
        Sig(:,:,k) = U*max(Sigma,0)*U';
    end
    
    mupo(:,:,n1)    = mu;
    Sigpo(:,:,:,n1) = Sig;
    wpo(:,:,n1)     = w;
    
    [~,ix] = sort(mu(1,:),2);

    mupo(:,:,n1)    = mu(:,ix);
    Sigpo(:,:,:,n1) = Sig(:,:,ix);
    wpo(1,:,n1)     = w(ix);    
end

normmupo = sqrt(sum(mupo.^2,1));
[~,mn]   = min(normmupo,[],3);
[~,mx]   = max(normmupo,[],3);

mn = mupo(:,:,mn(1));
mx = mupo(:,:,mx(end));

m0 = [];
for i=1:size(mn,1)
  m0 = [m0; linspace(mn(i,1),mx(i,end),K)];  
end

% n2 = 1; % Random subject for setting posteriors of all subjects

for n1=1:N
    cp(n1).pr.m = m0;
    cp(n1).pr.b = ones(1,K);
    for k=1:K
        cp(n1).pr.W(:,:,k) = eye(D,D);
    end
    cp(n1).pr.n = (D + 1)*ones(1,K);

    % Use 'responsibilities' from initialization to set sufficient statistics
    Nk   = Nf(n1)*wpo(:,:,n1);
    xbar = mupo(:,:,n1);
    S    = Sigpo(:,:,:,n1);
    
    b = cp(n1).pr.b + Nk;
    n = cp(n1).pr.n + Nk;
    W = zeros(D,D,K);
    for k = 1:K
        m(:,k) = (cp(n1).pr.b(k)*cp(n1).pr.m(:,k) + Nk(k).*xbar(:,k))./b(k);

        W0inv    = inv(cp(n1).pr.W(:,:,k));
        mlt1     = cp(n1).pr.b(k).*Nk(k)/(cp(n1).pr.b(k) + Nk(k));
        diff1    = xbar(:,k) - cp(n1).pr.m(:,k);
        W(:,:,k) = inv(W0inv + Nk(k)*S(:,:,k) + mlt1*(diff1*diff1'));
    end  

    cp(n1).po.m = m;
    cp(n1).po.b = b;
    cp(n1).po.W = W;
    cp(n1).po.n = n;

    cp(n1).dow = do.w;
    cp(n1).lnw = log(ones(1,K)/K);    
end

% for n=1:N  
%     %----------------------------------------------------------------------
%     % mixing weights
%     cp(n).dow = do.w;
%     cp(n).w   = ones(1,K);
%         
%     %----------------------------------------------------------------------
%     % Flat priors
%     
%     m0    = zeros(C,K);
%     beta0 = ones(1,K);   
%     nu0   = C*ones(1,K);
%     for k=1:K
%         W0(:,:,k) = eye(C,C);
%     end
% 
%     cp(n).pr.m    = m0;
%     cp(n).pr.beta = beta0;
%     cp(n).pr.W    = W0;                                                       
%     cp(n).pr.nu   = nu0;
%     
%     %----------------------------------------------------------------------
%     % Posteriors
%         
%     if isempty(pthmu)
%         m    = mupr;
%         beta = beta0;   
%         nu   = nu0;
%         for k=1:K
%             W(:,:,k) = nu0(k)*inv(Sigpr(:,:,k));
%         end  
%     else
%         %------------------------------------------------------------------
%         % Compute sufficient statistics from atlas
%         
%         f = zeros([prod(d) C]);
%         for c=1:C
%             Nii = nifti(pthf{n,c});
%             fc  = Nii.dat(:,:,:);
% 
%             fc(~cp(n).msk{c}) = 0;
% 
%             f(:,c) = fc(:);
%         end
%         
%         s0 = zeros(1,K);
%         s1 = zeros(C,K);
%         S2 = zeros(C,C,K);
%         for k=1:K
%             b = mu(:,:,:,k);
% 
%             s0(k)     = sum(sum(sum(b)));            
%             s1(:,k)   = sum(bsxfun(@times,f,b))/s0(k);
%             diff1     = bsxfun(@minus,f,s1(:,k)');
%             S2(:,:,k) = sum(b.*diff1.^2)/s0(k);
%         end   
% 
%         %------------------------------------------------------------------
%         % Compute posteriors from above sufficient statistics
% 
%         beta = beta0 + s0;
%         nu   = nu0 + s0;
%         m    = (repmat(beta0,C,1).*m0 + (ones(C,1)*s0).*s1)./(ones(C,1)*beta);
%         W    = zeros(C,C,K);
%         for k = 1:K
%             mult1    = beta0(k).*s0(k)/(beta0(k) + s0(k));
%             diff1    = s1(:,k) - m0(:,k);
%             W(:,:,k) = inv(inv(W0(:,:,k)) + s0(k)*S2(:,:,k) + mult1*(diff1*diff1'));
%         end 
%         
%         cp(n).pr.m    = m;
%         cp(n).pr.beta = beta;
%         cp(n).pr.W    = W;                                                       
%         cp(n).pr.nu   = nu;
%     end
%     
%     cp(n).po.m    = m;
%     cp(n).po.beta = beta;
%     cp(n).po.W    = W;
%     cp(n).po.nu   = nu;
% end
%==========================================================================