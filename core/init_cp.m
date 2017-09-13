function cp = init_cp(pthf,pthmu,mu,K,d,runpar)
[N,D] = size(pthf);

cp(N) = struct;
% parfor (n1=1:N,runpar)
for n1=1:N
    f = zeros([prod(d) D]);
    for c=1:D
        Nii    = nifti(pthf{n1,c});
        f(:,c) = reshape(Nii.dat(:,:,:),[],1);
    end
    Nf = nnz(f(:))/D;          

    % k-means--------------------------------------------------------------        
    f(f==0) = NaN;    
    
%     init{1} = 'eqspace';
%     init{2} = [];
%     labels  = Kmeans(f,K,init);   

    labels = label_data(f,K,d);
    
    % Get estimate of MoG parameters from suffstats------------------------
    mom        = SuffStats(f,labels);
    [w,mu,Sig] = GaussiansFromSuffStats(mom);

%     % To ensure positive semi-definite covariance matrix
%     for k=1:K
%         B          = (Sig(:,:,k) + Sig(:,:,k)')/2; 
%         [U,Sigma]  = eig(B); 
%         Sig(:,:,k) = U*max(Sigma,0)*U';
%     end
    
    [~,ix] = sort(mu(1,:),2);

    mu  = mu(:,ix);
    Sig = Sig(:,:,ix);
    w   = w(ix);    

    cp(n1).pr.m = zeros(D,K);
    cp(n1).pr.b = ones(1,K);
    for k=1:K
        cp(n1).pr.W(:,:,k) = 200*eye(D,D);
    end
    cp(n1).pr.n = 20*ones(1,K);

    % Use 'responsibilities' from initialization to set sufficient statistics
    s0 = Nf*w;
    s1 = mu;
    S2 = Sig;
    
    b = cp(n1).pr.b + s0;
    n = cp(n1).pr.n + s0;
    m = zeros(D,K);
    W = zeros(D,D,K);
    for k = 1:K
        m(:,k) = (cp(n1).pr.b(k)*cp(n1).pr.m(:,k) + s0(k).*s1(:,k))./b(k);

        W0inv    = inv(cp(n1).pr.W(:,:,k));
        mlt1     = cp(n1).pr.b(k).*s0(k)/(cp(n1).pr.b(k) + s0(k));
        diff1    = s1(:,k) - cp(n1).pr.m(:,k);
        W(:,:,k) = inv(W0inv + s0(k)*S2(:,:,k) + mlt1*(diff1*diff1'));
    end  

    % Init cluster struct
    cp(n1).po.m = m;
    cp(n1).po.b = b;
    cp(n1).po.W = W;
    cp(n1).po.n = n;
    cp(n1).dow  = 0;
    cp(n1).lnw  = log(ones([1,K],'single')/K);    
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

%==========================================================================
function mom = SuffStats(X,Q,mom)
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
    code = bitor(code,bitshift(feval(cast,isfinite(X(:,i)) & (X(:,i)~=0)),(i-1)));
end
for i=2:numel(mom),
    msk0      = mom(i).ind;
    ind       = find(code==msk0*(2.^(0:(M-1))'));
    if ~isempty(ind),
        x         = X(ind,msk0);
        for k=1:K,
            q                = Q(ind,k);
            mom(i).s0(1,k)   = mom(i).s0(1,k)   + sum(q);
            mom(i).s1(:,k)   = mom(i).s1(:,k)   + x'*q;
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + bsxfun(@times,q,x)'*x;
        end
    end
end
%==========================================================================

%==========================================================================
function [w,mn,vr] = GaussiansFromSuffStats(mom)
K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);
mg = zeros(1,K);
mn = zeros(N,K);
vr = zeros(N,N,K);
ll = -Inf;
for k=1:K,
    C  = eye(N);
    mu = zeros(N,1);
    s0 = 0;
    for i=1:numel(mom), s0 = s0 + mom(i).s0(k); end
    for iter=1:1024,
        s1     = zeros(N,1);
        S2     = zeros(N);
        P      = inv(C);
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_ll = ll;
        ll     = 0;
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

                % Compute objective function
                S              = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                ll             = ll - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi)))) - 0.5*trace(S/C(ind,ind));
            end
        end
        mu = s1/s0;
        C  = (S2 - s1*s1'/s0)/(s0+eps);
        %fprintf('%d\t%g\n', iter,ll);
        if ll-old_ll < 1e-12, break; end
    end
    mg(1,k)   = s0;
    mn(:,k)   = mu;
    vr(:,:,k) = C;
end
w = mg/sum(mg);
% fprintf('ML (simple)\t%g\n', ll);
% disp(mu)
% disp(C)
%==========================================================================

%==========================================================================
function labels = label_data(f,K,d)

opts = statset('MaxIter',1000);

labels = kmeans(f,K,'Options',opts);

labels                    = labels';
labels(~isfinite(labels)) = K + 1;

nlabels = zeros(prod(d),K + 1);

idx          = sub2ind(size(nlabels),1:prod(d),labels);
nlabels(idx) = 1;

idx = nlabels(:,K + 1) == 1;    
for k=1:K
   nlabels(idx,k) = NaN;
end
nlabels(:,K + 1) = [];

labels  = nlabels;    
%==========================================================================