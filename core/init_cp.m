function cp = init_cp(pthx,pthmsk,pthmu,mu,K,C)
N = numel(pthx);

init     = [];
init{1}  = 'rand';
init{2}  = [];

mupo  = zeros(C,K,N);
Sigpo = zeros(C,C,K,N);
for n=1:N
    Nii = nifti(pthx{n});
    f   = Nii.dat(:,:,:);

    Nii = nifti(pthmsk{n});
    msk = logical(Nii.dat(:,:,:));

    f = f(msk);

    [mupo(:,:,n),Sigpo(:,:,:,n)] = Kmeans(f,K,init);

    [mupo(:,:,n),ix] = sort(mupo(:,:,n),2);
    Sigpo(:,:,:,n)   = Sigpo(:,:,ix,n);
end

mupr  = mean(mupo,3);
Sigpr = mean(Sigpo,4);
        
for n=1:N  
    %----------------------------------------------------------------------
    % mixing weights
    cp(n).w = ones(1,K)/K;
        
    %----------------------------------------------------------------------
    % priors
    m0    = mupr;
    beta0 = 0.01*ones(1,K);   
    nu0   = C*ones(1,K) - 0.99;
    for k=1:K
        W0(:,:,k) = inv(Sigpr(:,:,k))/nu0(k);
    end
    
    cp(n).pr.m    = m0;
    cp(n).pr.beta = beta0;
    cp(n).pr.W    = W0;                                                       
    cp(n).pr.nu   = nu0;
    
    %----------------------------------------------------------------------
    % posteriors
        
    if isempty(pthmu)
        m    = mupo(:,:,n);
        beta = beta0;   
        nu   = nu0;
        for k=1:K
            W(:,:,k) = inv(Sigpo(:,:,k,n))/nu(k);
        end  
    else
        %------------------------------------------------------------------
        % Compute sufficient statistics from atlas

        Nii = nifti(pthmsk{n});
        msk = logical(Nii.dat(:,:,:));
        msk = msk(:);
        
        Nii = nifti(pthx{n});
        x   = Nii.dat(:,:,:);
        x   = x(msk);       
        
        s0 = zeros(1,K);
        s1 = zeros(C,K);
        S2 = zeros(C,C,K);
        for k=1:K
            b = mu(:,:,:,k); % Wrong
            b = b(msk);

            s0(k)     = sum(sum(sum(b)));            
            s1(:,k)   = sum(bsxfun(@times,x,b))/s0(k);
            diff1     = bsxfun(@minus,x,s1(:,k)');
            S2(:,:,k) = sum(b.*diff1.^2)/s0(k);
        end   

        %------------------------------------------------------------------
        % Compute posteriors from above sufficient statistics

        beta = beta0 + s0;
        nu   = nu0 + s0;
        m    = (repmat(beta0,C,1).*m0 + (ones(C,1)*s0).*s1)./(ones(C,1)*beta);
        W    = zeros(C,C,K);
        for k = 1:K
            mult1    = beta0(k).*s0(k)/(beta0(k) + s0(k));
            diff1    = s1(:,k) - m0(:,k);
            W(:,:,k) = inv(inv(W0(:,:,k)) + s0(k)*S2(:,:,k) + mult1*(diff1*diff1'));
        end 
        
        cp(n).pr.m    = m;
        cp(n).pr.beta = beta;
        cp(n).pr.W    = W;                                                       
        cp(n).pr.nu   = nu;
    end
    
    cp(n).po.m    = m;
    cp(n).po.beta = beta;
    cp(n).po.W    = W;
    cp(n).po.nu   = nu;
end
%==========================================================================

%==========================================================================
function [mg,mu,sig] = fit_gmm2h(h,x,K,shwinfo)
if nargin < 4, shwinfo = 0; end

mg  = ones(K,1)/K;
mu  = linspace(min(x),max(x),K)'./K;
sig = ones(K,1)*(max(x) - min(x))./K;  

m0    = zeros(K,1);
m1    = zeros(K,1);
m2    = zeros(K,1);
ll(1) = -Inf;
for iter=1:10000,
p  = zeros(numel(x),K);
for k=1:K,
    % Product Rule
    % p(class=k, intensity | mg, nu, sig) = p(class=k|mg) p(intensity | nu, sig, class=k)
    p(:,k) = mg(k)*normpdf(x(:),mu(k),sig(k));
end

% Sum Rule
% p(intensity | mg, nu, sig) = \sum_k p(class=k, intensity | mg, nu, sig)
sp         = sum(p,2)+eps;
ll(iter+1) = sum(log(sp).*h(:));
if ll(iter+1) - ll(iter) < 1e-8*sum(h)
    if shwinfo == 2,
        figure(4001);
        md = mean(diff(x));
        plot(x(:),(h/sum(h))/md,'b-',x(:),sp,'r-'); hold on
        plot(x(:),p,'--');        
        set(gca,'TickLabelInterpreter','latex');  
        xlabel('Image intensity','Interpreter','latex')
        ylabel('Probability','Interpreter','latex')
        legend({'Empirical','Fit','Air','Tissue'},'Interpreter','latex');
        drawnow;
    end
    break; 
end

if shwinfo == 3,
    figure(4001);
    subplot(121); plot(0:numel(ll)-2,ll(2:end))  
    md = mean(diff(x));
    subplot(122); plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); 
    drawnow
end

% Bayes Rule
% p(class=k | intensity, mg, nu, sig) = p(class=k, intensity | mg, nu, sig) / p(intensity | mg, nu, sig)
p = bsxfun(@rdivide,p,sp);

% Compute moments from the histograms, weighted by the responsibilities (p).
for k=1:K,
    m0(k) = sum(p(:,k).*h(:));             % Number of voxels in class k
    m1(k) = sum(p(:,k).*h(:).*x(:));       % Sum of the intensities in class k
    m2(k) = sum(p(:,k).*h(:).*x(:).*x(:)); % Sum of squares of intensities in class k
end
mg = m0/sum(m0);
for k=1:K,
    mu(k) = m1(k)./m0(k);                                 % Mean
    sig(k) = (m2(k)-m1(k)*m1(k)/m0(k)+1e-6)/(m0(k)+1e-6); % Variance
end
sig = sqrt(sig);
end
%==========================================================================   