function obj = spm_segment(obj,pth_template,fig)

% Load template
%-----------------------------------------------------------------------
tpm = spm_load_priors8(pth_template);

% Parameters, etc.
%-----------------------------------------------------------------------
Affine    = obj.Affine;
V         = obj.image;
N         = numel(V);
d0        = V(1).dim(1:3);
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(obj.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
d         = [size(x0) length(z0)];
lkp       = obj.lkp;
Kb        = max(obj.lkp);
K_lab     = obj.K_lab;
tol1      = obj.tol1;
wp_reg    = 1;
modality  = obj.modality;
niter     = obj.niter;
nitgmm    = obj.nitgmm;
nsubit    = obj.nsubit;
do_bf     = obj.do_bf;
if obj.uniform, obj.do_def = false; end
do_def    = obj.do_def;
do_wp     = obj.do_wp;
print_ll  = obj.print_ll;
pth_vel   = obj.pth_vel;
M         = tpm.M\Affine*V(1).mat;
tot_S     = obj.tot_S;

% Some random numbers are used, so initialise random number generators to
% give the same results each time.
%-----------------------------------------------------------------------
rng('default');
rng(1);

% Fudge Factor - to (approximately) account for non-independence of voxels.
%-----------------------------------------------------------------------
ff = compute_fudge_factor(obj,V,sk);

% Load data into buffer
%-----------------------------------------------------------------------
[buf,nm,vr0,mn,mx] = init_buf(N,obj,V,x0,y0,z0,o,M,tpm,tot_S);

% Initialise weights
%-----------------------------------------------------------------------
if isfield(obj,'wp'), wp = obj.wp;
else,                 wp = ones(1,Kb)/Kb;
end

wp(K_lab{2}) = eps;
wp           = wp/sum(wp);
        
if isfield(obj,'mg')
    mg  = obj.mg;     
else              
    mg  = ones(Kb,1);
    lkp = 1:Kb;   
end

% Initialise bias field
%-----------------------------------------------------------------------
[buf,chan,llrb] = init_bf(buf,N,obj,V,x0,y0,z0,ff);

% Initialise deformation and template
%-----------------------------------------------------------------------
[buf,param,MT,sk4,Twarp,llr] = init_def(buf,obj,lkp,sk,vx,ff,d,fig,wp,x0,y0,z0,tpm,M);

% Initialise GMM
%-----------------------------------------------------------------------
gmm = init_gmm(obj,N,buf,vr0,mn,mx);                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iterating
%------------------------------------------------------------
armijo = ones(1,2);
ll     = -Inf;
L      = {ll,llrb,llr};
for iter=1:niter

    for subit=1:nsubit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate cluster parameters
        %------------------------------------------------------------
        [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,K_lab,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll);
 
        debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,K_lab);
        
        if subit > 1 && ~((ll-ooll)>2*tol1*nm), break; end
        ooll = ll;

        if do_bf
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate bias
            % Note that for multi-spectral data, the covariances among
            % channels are not computed as part of the second derivatives.
            % The aim is to save memory, and maybe make the computations
            % faster.
            %------------------------------------------------------------
            [ll,llrb,buf,chan,L,armijo(1)] = update_bf(ll,llrb,llr,buf,mg,gmm,wp,lkp,K_lab,chan,fig,L,print_ll,armijo(1));
        
            debug_view('bf',fig{2},lkp,buf,modality);
        end
        
        if iter==1 && subit==1
            % Most of the log-likelihood improvements are in the first iteration.
            % Show only improvements after this, as they are more clearly visible.
            spm_plot_convergence('Clear');
            spm_plot_convergence('Init','Processing','Log-likelihood','Iteration');

            if numel(obj.lkp) ~= numel(lkp)
                lkp      = obj.lkp;
                K        = numel(lkp);
                Kb       = max(lkp);
                [gmm,mg] = more_gmms(gmm,lkp,K,Kb);
            end
        end
    end

    if do_def
        for subit=1:nsubit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate cluster parameters
            %------------------------------------------------------------
            [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,K_lab,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll);

            debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,K_lab);

            if subit > 1 && ~((ll-oll)>2*tol1*nm), break; end
            oll = ll;
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate deformations
            %------------------------------------------------------------
            [ll,llr,buf,Twarp,L,armijo(2)] = update_def(ll,llrb,llr,buf,mg,gmm,wp,lkp,K_lab,Twarp,sk4,M,MT,tpm,x0,y0,z0,param,iter,fig,L,print_ll,tot_S,armijo(2));
        end    
    end

    debug_view('template',fig{3},lkp,buf,wp);  
    
    if iter>=10 && ~((ll-ooll)>2*tol1*nm)
        % Finished
        break
    end
end
clear buf tpm

% For setting the DC component of all the bias fields so that they
% average to 0 (used for global signal normalisation).
%--------------------------------------------------------------------------
bf_dc               = zeros(1,N);
for n=1:N, bf_dc(n) = chan(n).T(1,1,1); end 

% Save the results
%--------------------------------------------------------------------------
Nii = nifti(pth_vel);
for i=1:3
    Nii.dat(:,:,:,i) = Twarp(:,:,:,i); 
end
clear Twarp Nii

obj.Affine = Affine;
obj.lkp    = lkp;
obj.MT     = MT;
obj.Tbias  = {chan(:).T};
obj.wp     = wp;
obj.mg     = mg;
obj.gmm    = gmm;
obj.ll     = ll;
obj.new    = true;
obj.bf_dc  = bf_dc;
obj.d0     = d0;
return;
%=======================================================================

%=======================================================================
function ff = compute_fudge_factor(obj,V,sk)
% Fudge Factor - to (approximately) account for non-independence of voxels.
% Note that variances add, and that Var[a*x + b*y] = a^2*Var[x] + b^2*Var[y]
% Therefore the variance of i.i.d. noise after Gaussian smoothing is equal
% to the sum of the Gaussian function squared times the original variance.
% A Gaussian is given by g=sqrt(2*pi*s^2)^(-1/2)*exp(-0.5*x.^2/s^2);
% After squaring, this is (2*pi*s^2)^(-1)*exp(-x.^2/s^2), which is a scaled
% Gaussian. Letting s2 = 2/sqrt(2), this is equal to
% (4*pi*s^2)^(-1/2)*(2*pi*s2^2)^(-1/2)*exp(-0.5*x.^2/s2^2), from which
% the (4*pi*s^2)^(-1/2) factor comes from.
fwhm = obj.fwhm;                            % FWHM of image smoothness
vx   = sqrt(sum(V(1).mat(1:3,1:3).^2));     % Voxel size
fwhm = fwhm+mean(vx); 
s    = fwhm/sqrt(8*log(2));                 % Standard deviation
ff   = prod(4*pi*(s./vx./sk).^2 + 1)^(1/2); 
%=======================================================================

%=======================================================================
function [gmm,mg] = more_gmms(gmm,lkp,K,Kb)
mg = ones(K,1)/K;

if gmm.ml
    mn1 = gmm.mn;
    vr1 = gmm.vr;
    N   = size(mn1,1);

    % Use moments to compute means and variances, and then use these
    % to initialise the Gaussians
    mn = ones(N,K);
    vr = zeros(N,N,K);

    for k=1:Kb
        % A crude heuristic to replace a single Gaussian by a bunch of Gaussians
        % If there is only one Gaussian, then it should be the same as the
        % original distribution.
        kk  = sum(lkp==k);
        w   = 1./(1+exp(-(kk-1)*0.25))-0.5;
        
        mn(:,lkp==k)   = sqrtm(vr1(:,:,k))*randn(N,kk)*w + repmat(mn1(:,k),[1,kk]);
        vr(:,:,lkp==k) = repmat(vr1(:,:,k)*(1-w),[1,1,kk]);
        
        mg(lkp==k)     = 1/kk;
    end

    gmm.mn = mn;
    gmm.vr = vr;
else
    % Adjust priors
    %----------------------------------------------------------------------
    m0 = gmm.pr.m;
    n0 = gmm.pr.n;
    b0 = gmm.pr.b;
    W0 = gmm.pr.W;
    
    m = zeros(size(m0));
    n = zeros(size(n0));
    b = zeros(size(b0));
    W = zeros(size(W0));
    
    for k=1:Kb
        kk            = sum(lkp==k);
        b(lkp==k)     = b0(k);
        m(:,lkp==k)   = repmat(m0(:,k),[1,kk]);
        n(lkp==k)     = n0(k);
        W(:,:,lkp==k) = repmat(W0(:,:,k),[1 1 kk]); 
        
        mg(lkp==k) = 1/kk;
    end
    
    gmm.pr.m = m;
    gmm.pr.n = n;
    gmm.pr.b = b;
    gmm.pr.W = W;
    
    % Adjust posteriors
    %----------------------------------------------------------------------    
    m0 = gmm.po.m;
    n0 = gmm.po.n;
    b0 = gmm.po.b;
    W0 = gmm.po.W;
 
    m = zeros(size(m0));
    n = zeros(size(n0));
    b = zeros(size(b0));
    W = zeros(size(W0));

    N = size(m,1);
    
    for k=1:Kb
        kk  = sum(lkp==k);
        w   = 1./(1+exp(-(kk-1)*0.25))-0.5;
        
        vr0           = inv(n0(k)*W0(:,:,k));
        pr0           = inv(vr0*(1 - w));                        
        b(lkp==k)     = b0(k);
        m(:,lkp==k)   = sqrtm(vr0)*randn(N,kk)*w + repmat(m0(:,k),[1,kk]);
        n(lkp==k)     = n0(k);
        W(:,:,lkp==k) = repmat(pr0/n0(k),[1 1 kk]);   
    end
    
    gmm.po.m = m;
    gmm.po.n = n;
    gmm.po.b = b;
    gmm.po.W = W;
end
%=======================================================================  