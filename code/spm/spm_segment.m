function obj = spm_segment(obj,fig)

% Load template
%-----------------------------------------------------------------------
tpm = spm_load_logpriors(obj.pth_template);

% Parameters, etc.
%-----------------------------------------------------------------------
V         = obj.image; 
N         = numel(V);
d0        = V(1).dim(1:3);
vx        = vxsize(V(1).mat);
sk        = max([1 1 1],round(obj.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
d         = [size(x0) length(z0)];
lkp       = obj.lkp;
Kb        = max(lkp.part);
tol1      = obj.tol1;
wp_reg    = 1;
wp_lab    = obj.wp_lab;
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
Affine    = obj.Affine;
M         = tpm.M\Affine*V(1).mat;
tot_S     = obj.tot_S;

% Fudge Factor - to (approximately) account for non-independence of voxels.
%-----------------------------------------------------------------------
ff = compute_fudge_factor(obj.fwhm,vx,sk);

% Load data into buffer
%-----------------------------------------------------------------------
[buf,nm,vr0,mn,mx,scl_int] = init_buf(N,obj,V,x0,y0,z0,o,M,tpm,tot_S);

% Initialise weights
%-----------------------------------------------------------------------
if isfield(obj,'wp'), wp = obj.wp;
else,                 wp = ones(1,Kb)/Kb;
end
        
part0 = lkp.part;
if isfield(obj,'mg')
    mg = obj.mg;   
    if numel(mg)~=numel(lkp.part)
        lkp.part = 1:Kb;       
    end
else              
    mg       = ones(Kb,1);
    lkp.part = 1:Kb;   
end

% Initialise bias field
%-----------------------------------------------------------------------
[buf,chan,llrb] = init_bf(buf,N,obj,V,x0,y0,z0,ff,scl_int);

% Initialise deformation and template
%-----------------------------------------------------------------------
[buf,param,MT,sk4,Twarp,llr] = init_def(buf,obj,sk,vx,ff,d,fig,wp,x0,y0,z0,tpm,M);

% Initialise GMM
%-----------------------------------------------------------------------
[gmm,buf] = init_gmm(obj,N,buf,vr0,mn,mx);                               

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
        [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_lab);
 
        debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,wp_lab);
        
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
            [ll,llrb,buf,chan,L,armijo(1)] = update_bf(ll,llrb,llr,buf,mg,gmm,wp,lkp,chan,fig,L,print_ll,armijo(1),wp_lab);
        
            debug_view('bf',fig{2},lkp,buf,modality);
        end
        
        if iter==1 && subit==1
            % Most of the log-likelihood improvements are in the first iteration.
            % Show only improvements after this, as they are more clearly visible.
            spm_plot_convergence('Clear');
            spm_plot_convergence('Init','Processing','Log-likelihood','Iteration');

            if numel(lkp.part) ~= numel(part0)
                lkp.part = part0;
                K        = numel(lkp.part);
                Kb       = max(lkp.part);
                [gmm,mg] = more_gmms(gmm,lkp.part,K,Kb);
            end
        end
    end

    if do_def
        for subit=1:nsubit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate cluster parameters
            %------------------------------------------------------------
            [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_lab);

            debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,wp_lab);

            if subit > 1 && ~((ll-oll)>2*tol1*nm), break; end
            oll = ll;
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate deformations
            %------------------------------------------------------------
            [ll,llr,buf,Twarp,L,armijo(2)] = update_def(ll,llrb,llr,buf,mg,gmm,wp,lkp,Twarp,sk4,M,MT,tpm,x0,y0,z0,param,iter,fig,L,print_ll,tot_S,armijo(2),wp_lab);
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
function [gmm,mg] = more_gmms(gmm,part,K,Kb)
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
        kk  = sum(part==k);
        w   = 1./(1+exp(-(kk-1)*0.25))-0.5;
        
        mn(:,part==k)   = sqrtm(vr1(:,:,k))*randn(N,kk)*w + repmat(mn1(:,k),[1,kk]);
        vr(:,:,part==k) = repmat(vr1(:,:,k)*(1-w),[1,1,kk]);
        
        mg(part==k)     = 1/kk;
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
        kk            = sum(part==k);
        b(part==k)     = b0(k);
        m(:,part==k)   = repmat(m0(:,k),[1,kk]);
        n(part==k)     = n0(k);
        W(:,:,part==k) = repmat(W0(:,:,k),[1 1 kk]); 
        
        mg(part==k) = 1/kk;
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
        kk  = sum(part==k);
        w   = 1./(1+exp(-(kk-1)*0.25))-0.5;
        
        vr0           = inv(n0(k)*W0(:,:,k));
        pr0           = inv(vr0*(1 - w));                        
        b(part==k)     = b0(k)/kk;
        m(:,part==k)   = sqrtm(vr0)*randn(N,kk)*w + repmat(m0(:,k),[1,kk]);
        n(part==k)     = n0(k)/kk;
        W(:,:,part==k) = repmat(pr0/n0(k),[1 1 kk]);   
    end
    
    gmm.po.m = m;
    gmm.po.n = n;
    gmm.po.b = b;
    gmm.po.W = W;
end
%=======================================================================  