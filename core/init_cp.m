function cp = init_cp(Px,Pmsk,Pmu,mu,K,C,tempdir,mat,d)
N = numel(Px);

mu  = zeros(C,K,N);
Sig = zeros(C,C,K,N);
if isempty(Pmu)
    for n=1:N
        Nii = nifti(Px{n});
        f   = Nii.dat(:,:,:);

        Nii = nifti(Pmsk{n});
        msk = logical(Nii.dat(:,:,:));

        f = f(msk);

        init     = [];
        init{1}  = 'rand';
        init{2}  = [];
        [mu(:,:,n),Sig(:,:,:,n)] = Kmeans(f,K,init);
    end
    
    mu  = mean(mu,3);
    Sig = mean(Sig,4);
end    
        
for n=1:N
    Nii = nifti(Pmsk{n});
    msk = logical(Nii.dat(:,:,:));

    Nii = nifti(Px{n});
    x   = single(Nii.dat(:,:,:));
    x   = x(msk);

    %------------------------------------------------------------------
    % Set priors
    m0    = mu;
    beta0 = ones(1,K);   
    nu0   = 20*ones(1,K);
    for k=1:K
        W0(:,:,k) = inv(Sig(:,:,k))/nu0(k);
    end
    
    %----------------------------------------------------------------------
    % Set posteriors
        
    if isempty(Pmu)
        m    = m0;
        beta = beta0;        
        nu   = nu0;
        W    = W0;
    else
        %------------------------------------------------------------------
        % Compute sufficient statistics from atlas

        s0 = zeros(1,K);
        s1 = zeros(C,K);
        S2 = zeros(C,C,K);
        for k=1:K
            b = mu(:,:,:,k);
            b = b(msk);

            s0(1,k)   = sum(sum(sum(b)));            
            s1(:,k)   = sum(bsxfun(@times,x,b))/s0(1,k);
            diff      = bsxfun(@minus,x,s1(:,k)');
            S2(:,:,k) = 1/s0(k).*sum(b.*diff.^2);
        end   

        %------------------------------------------------------------------
        % Compute posteriors from above sufficient statistics

        beta = beta0 + s0;
        nu   = nu0 + s0;
        m    = (repmat(beta0,C,1).*m0 + (ones(C,1)*s0).*s1)./(ones(C,1)*beta);
        W    = zeros(C,C,K);
        for k = 1:K
            mult1    = beta0(k).*s0(k)/(beta0(k) + s0(k));
            diff3    = s1(:,k) - m0(:,k);
            W(:,:,k) = inv(inv(W0(:,:,k)) + s0(k)*S2(:,:,k) + mult1*(diff3*diff3'));
        end 
    end
    
    %------------------------------------------------------------------
    % Construct cluster struct

    cp(n).w = ones(1,K)/K;

    cp(n).pr.m    = m0;
    cp(n).pr.beta = beta0;
    cp(n).pr.W    = W0;                                                       
    cp(n).pr.nu   = nu0;

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

%==========================================================================
function [varargout] = Kmeans(xv,nC,init,varargin)   
xv(xv==0) = NaN;

label = NaN(size(xv,1),nC);

N  = size(xv,1);
D  = size(xv,2);
M  = 2^D;   

in = bi2de(double(isfinite(xv)))+1;
%initialize means------------------------------------------------------

if strcmp(init{1},'eqspace')
    if isempty(init{2})
    st = [1/nC:1/nC:1]'; 
    else
        o  = init{2};
        st = [o:(1-o)/(nC-1):1]'; 
    end
    mu = permute(0.8*st*max(xv),[3,2,1]);

elseif strcmp(init{1},'rand')
     if ~isempty(init{2})
         xv(sum(xv,2)<=init{2}) = NaN;
         in = bi2de(double(isfinite(xv))) + 1;
     end
    mu  = permute(xv(randsample(find(in==M),nC),:),[3,2,1]);

elseif strcmp(init{1},'custom')
    if isempty(init{2})
        error('Have to specify means initialization')
    else
        st = init{2};
        mu = permute(st'*max(xv),[3,2,1]);
    end
end

for m=2:M

    conf_log = logical(de2bi(m-1,D));

    xvsub  = xv(repmat(m,[N,1])==in,conf_log);
    musub  = mu(1,conf_log,:);        
    N1     = size(xvsub,1); 
    D1     = size(xvsub,2);        
    dis    = zeros(N1,nC);

    %start iterations--------------------------------------------------               
    for it=1:20

        %labeling------------------------------------------------------
        for k=1:nC
        dis(:,k) = squeeze(sqrt(sum(bsxfun(@minus,xvsub,musub(:,:,k)).^2,2))); 
        end
        labelsub = bsxfun(@minus,dis,min(dis,[],2))==0;

        %updating means------------------------------------------------
        labelw = permute(repmat(labelsub,[1,1,D1]),[1,3,2]);
        m0     = sum(labelw,1);
        m1     = sum(labelw.*repmat(xvsub,[1,1,nC]),1);
        clear labelw
        musub  = m1./m0;
    end
    clear dis

    label(repmat(m,[N,1])==in,:) = labelsub;
    clear labelsub
end
clear xvsub

mom = spm_SuffStats(xv,label);
if isempty(varargin)
    [~,varargout{1},varargout{2},~] = spm_VBGaussiansFromSuffStats_v2(mom);
else
    varargout{1} = spm_VBGaussiansFromSuffStats_v2(mom,varargin{1});
end
%==========================================================================                 