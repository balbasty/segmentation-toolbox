function [obj,munum,muden] = spm_preprocx(obj,logtpm)
% New and improved (?) unified segmentation routine
%
% FORMAT [obj,munum,muden] = spm_preprocx(obj,logtpm)
%
%_______________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

V         = obj.image;
descrip   = obj.descrip;
N         = numel(V);
Affine    = obj.Affine;
M         = logtpm.M\Affine*V(1).mat;
d0        = V(1).dim(1:3);
dtpm      = logtpm.d;
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(obj.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
wp_reg    = obj.wp_reg;
tol1      = obj.tolseg;
tiny      = eps*eps;
lkp       = obj.lkp;
if isempty(lkp)
    K       = 2000;
    Kb      = numel(logtpm.dat);
    use_mog = false;
    vb      = false;
else
    K       = numel(obj.lkp);
    Kb      = max(obj.lkp);
    use_mog = true;
    vb      = obj.vb;
end
% lkpsemi = obj.lkpsemi;

kron = @(a,b) spm_krutil(a,b);

niter      = obj.niter;
niter1     = obj.niter1;
nsubitmog  = obj.nsubitmog;
nsubitbf   = obj.nsubitbf;
nitdef     = obj.nitdef;

dobias = obj.dobias;
dodef  = obj.dodef;
dotpm  = obj.dotpm;

% TODO: remove once missing data part works...
nomiss = 1;

% Some random numbers are used, so initialise random number generators to
% give the same results each time.
rng('default');
rng(1);

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

spm_diffeo('boundary',1);

% Initialise Deformation
%-----------------------------------------------------------------------
% This part is fiddly because of the regularisation of the warps.
% The fact that displacement fields are only parameterised every few
% voxels means that the functions in spm_diffeo need tweaking to
% account for the difference between the units of displacement and
% the separation of the voxels (if that makes sense).

% More work/thought is needed in terms of adjusting regularisation to
% account for different voxel sizes.  I'm still not satisfied that
% this (rescaling the regularisaiton by prod(vx.*sk)) is optimal.
% The same thing applies to all the nonlinear warping code in SPM.
param  = [sk.*vx prod(vx.*sk)*ff*obj.reg]; % FIX THIS (remove "prod(vx.*sk)")

% Mapping from indices of subsampled voxels to indices of voxels in image(s).
MT     = [sk(1) 0 0 (1-sk(1));0 sk(2) 0 (1-sk(2)); 0 0 sk(3) (1-sk(3));0 0 0 1];

% For multiplying and dividing displacements to map from the subsampled voxel indices
% and the actual image voxel indices.
sk4    = reshape(sk,[1 1 1 3]);

% Dimensions of sub-sampled image
d      = [size(x0) length(z0)];

if isfield(obj,'Twarp')
    Twarp = obj.Twarp; 
else
    Twarp = zeros([d,3],'single');
end
llr = -0.5*sum(sum(sum(sum(Twarp.*bsxfun(@times,spm_diffeo('vel2mom',bsxfun(@times,Twarp,1./sk4),param),1./sk4)))));

% Initialise bias correction
%-----------------------------------------------------------------------
cl   = cell(N,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
if use_mog
    chan = struct(args{:});
else
    chan = struct(args{:},'hist',cl,'lik',cl,'alph',cl,'grad',cl,'lam',cl,'interscal',cl);
end

for n=1:N
    % GAUSSIAN REGULARISATION for bias correction
    fwhm    = obj.biasfwhm(n);
    biasreg = obj.biasreg(n);
    vx      = sqrt(sum(V(n).mat(1:3,1:3).^2));
    d0      = V(n).dim;
    sd      = vx(1)*d0(1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
    sd      = vx(2)*d0(2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
    sd      = vx(3)*d0(3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
    Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*biasreg*ff;
    chan(n).C   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));

    % Basis functions for bias correction
    chan(n).B3  = spm_dctmtx(d0(3),d3(3),z0);
    chan(n).B2  = spm_dctmtx(d0(2),d3(2),y0(1,:)');
    chan(n).B1  = spm_dctmtx(d0(1),d3(1),x0(:,1));

    % Initial parameterisation of bias field
    if isfield(obj,'Tbias') && ~isempty(obj.Tbias{n})
        chan(n).T = obj.Tbias{n};
    else
        chan(n).T = zeros(d3);
    end
end

if isfield(obj,'msk') && ~isempty(obj.msk)
    VM = spm_vol(obj.msk);
    if sum(sum((VM.mat-V(1).mat).^2)) > 1e-6 || any(VM.dim(1:3) ~= V(1).dim(1:3))
        error('Mask must have the same dimensions and orientation as the image.');
    end
end

% Load the data
%-----------------------------------------------------------------------
im = zeros(N,1); % Intensity of voxels (for each channel)

scrand = zeros(N,1);
for n=1:N
    if spm_type(V(n).dt(1),'intt')
        scrand(n) = V(n).pinfo(1);
    end
end

% Overall moments used later for regularising via a ``Wishart-style prior''
s0 = zeros(1,N);
s1 = zeros(1,N);
S2 = zeros(1,N);

cl  = cell(length(z0),1);
buf = struct('f',cl,'mu',cl,'bf',cl,'code',cl);
for z=1:length(z0)        
    buf(z).msk = cell(1,N);
    if ~dotpm && d(3)>1
        % Load only those voxels that are more than 5mm up
        % from the bottom of the tissue probability map.  This
        % assumes that the affine transformation is pretty close.
	 
        %x1  = M(1,1)*x0 + M(1,2)*y0 + (M(1,3)*z0(z) + M(1,4));
        %y1  = M(2,1)*x0 + M(2,2)*y0 + (M(2,3)*z0(z) + M(2,4));
        z1  = M(3,1)*x0 + M(3,2)*y0 + (M(3,3)*z0(z) + M(3,4));
        e   = sqrt(sum(logtpm.M(1:3,1:3).^2));
        e   = 5./e; % mm from edge of TPM
        for n=1:N
            buf(z).msk{n} = z1>e(3);        
        end 
    else
        for n=1:N
            buf(z).msk{n} = ones(d(1:2),'logical');        
        end  
    end
        
    % Initially load all the data, but prepare to exclude
    % locations where any of the images is not finite, or
    % is zero.  We want this to work for skull-stripped
    % images too.
    fz = cell(1,N);
    for n=1:N
        fz{n}         = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        if nomiss
            buf(z).msk{1} = buf(z).msk{1} & get_msk(fz{n},descrip);            
        else
            buf(z).msk{n} = get_msk(fz{n},descrip);
        end
    end
    if nomiss
        for n=2:N
            buf(z).msk{n} = buf(z).msk{1};
        end
    end
    
    if isfield(obj,'msk') && ~isempty(obj.msk)
        % Exclude any voxels to be masked out
        msk        = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        for n=1:N
            buf(z).msk{n} = buf(z).msk{n} & msk;
        end
    end
    
    % Count voxels
    tmp        = 0;
    buf(z).nmn = cell(1,N);
    for n=1:N
        buf(z).nmn{n} = nnz(buf(z).msk{n});  % number of voxels in each channel
        tmp           = tmp + buf(z).msk{n}(:);
    end
    buf(z).nm = sum(tmp>0); % number of voxels in all channels    

    % Sum up intensities of all voxels (used for later normalising image
    % intensities)
    for n=1:N
        im(n) = im(n) + sum(fz{n}(buf(z).msk{n}));
    end
end
clear tmp

% For dealing with missing data
if N<=8
    cast = @uint8;
    typ  = 'uint8';
elseif N<=16
    cast = @uint16;
    typ  = 'uint16';
elseif N<=32
    cast = @uint32;
    typ  = 'uint32';
elseif N<=64
    cast = @uint64;
    typ  = 'uint64';
else
    error('Too many dimensions.');
end

for z=1:length(z0)  
    code = zeros([prod(d(1:2)) 1],typ);
    for n=1:N
        code = bitor(code,bitshift(feval(cast,buf(z).msk{n}(:)),(n-1)));
    end
    buf(z).code = code;
end
clear code

nm = 0;
for z=1:length(z0)
    buf(z).nm = nnz(buf(z).code);  % number of voxels per slice (over all channels)
    nm        = nm + buf(z).nm;    % number of voxels in all slices (over all channels)
end

a = zeros(1,N);
for z=1:length(z0)  
    for n=1:N                   
        fz = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        
        % Normalise image intensities (not if CT)            
        if strcmp(descrip,'CT')
            a(n) = 1.0;
        else
            a(n) = 512/(im(n)/nm);
        end
       
        % Eliminate unwanted voxels
        if scrand(n)
            % The variances used in a GMM may not be ideal because of the
            % discretisation of image intensities (Data is an integer
            % type), so to prevent aliasing in the histogram, small random
            % values are added.  It's not elegant, but the alternative
            % would be too slow for practical use.
            buf(z).f{n}  = single(a(n)*fz(buf(z).msk{n})+rand(buf(z).nmn{n},1)*scrand(n)-scrand(n)/2);
        else
            buf(z).f{n}  = single(a(n)*fz(buf(z).msk{n}));
        end    
        
        % Create moments to be used for constructing a variance prior
        s0(n) = s0(n) + buf(z).nmn{n};
        s1(n) = s1(n) + sum(buf(z).f{n});
        S2(n) = S2(n) + sum(buf(z).f{n}.^2);
    end

    % Create a buffer for tissue probability info
    buf(z).mu = zeros([buf(z).nm Kb],'single');
end

% Construct a ``Wishart-style prior''
% Mostly here to prevent some of the numerical instabilities that arise
% from computing the variance from the moments in a single pass.
vr0 = diag(S2./s0 - (s1./s0).^2)/Kb^2;
%for n=1:N
%    if spm_type(V(n).dt(1),'intt')
%        vr0(n,n) = vr0(n,n) + 0.083*V(n).pinfo(1,1);
%    end
%end

% Create initial bias field
%-----------------------------------------------------------------------
llrb = 0;
for n=1:N
    B1 = chan(n).B1;
    B2 = chan(n).B2;
    B3 = chan(n).B3;
    C  = chan(n).C;
    T  = chan(n).T;
    chan(n).ll = double(-0.5*T(:)'*C*T(:));
    for z=1:numel(z0)
        bf           = transf(B1,B2,B3(z,:),T);
        tmp          = bf(buf(z).msk{n});
        chan(n).ll   = chan(n).ll + double(sum(tmp));
        buf(z).bf{n} = single(exp(tmp));
    end
    llrb = llrb + chan(n).ll;
    clear B1 B2 B3 T C tmp
end

if isfield(obj,'wp')
    wp = obj.wp;
else
    wp = ones(1,Kb)/Kb;
end

% % Set weigths to zero and re-normalise (for semi-supervised approach)
% wp(lkpsemi) = eps;
% wp          = wp/sum(wp);

% For debugging
%-----------------------------------------------------------------------
zix     = floor(d(3)/2 + 1);
fig     = obj.fig;
verbose = obj.verbose;
slice   = zeros(d(1:2),'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run algorithm
%------------------------------------------------------------

spm_plot_convergence('Init','Initialising','Log-likelihood','Iteration');

ll = -Inf;
for iter=1:niter

    % Load the warped prior probability images into the buffer
    %------------------------------------------------------------
    for z=1:length(z0)
        if ~buf(z).nm, continue; end
        
        [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,buf(z).msk);        
        mu         = spm_sample_logpriors8(logtpm,x1,y1,z1);
        for k1=1:Kb
            buf(z).mu(:,k1) = mu{k1};
        end         
        clear x1 y1 z1
                
        if ~isempty(fig{2}) && z==zix && iter==1
            % Visualise initial template
            K1 = floor(sqrt(Kb));
            K2 = ceil(Kb/K1); 
            set(0,'CurrentFigure',fig{2});                                        
            for i=1:Kb            
                subplot(K1,K2,i);
                slice(buf(z).msk{1}) = mu{i};
                imagesc(slice'); axis image xy off; colormap(gray);
            end 
            drawnow
        end
        clear mu
    end
           
    % Starting estimates for intensity distribution parameters
    %------------------------------------------------------------
    if iter==1                  
        if use_mog
            % Starting estimates for Gaussian parameters
            %--------------------------------------------------------------
            mn = 0; vr = 0; po = 0;
                                
            K   = Kb;
            lkp = 1:Kb;                                                            
            mg  = ones(1,K);            
            if ~logtpm.uniform
                % Use data and template to estimate moments
                %------------------------------------------------------    
                mom = mom_struct(K,N,nomiss);  
                for z=1:length(z0)
                    mom = suffstats(mom,buf(z),buf(z).mu);
                end
            else
                % Use a k-means algorithm to estimate moments
                %------------------------------------------------------
                mom = spm_kmeans2mom(buf,K,nomiss,vr0);
            end

            if vb                                                
                % Use moments to compute initial Gaussian-Wishart
                % posteriors and priors
                %--------------------------------------------------
                [lq,~,po,pr] = spm_GaussiansFromSuffStats(vb,mom);
                
                if logtpm.uniform
                    % Sort according to m values of first channel
                    [~,ix] = sort(po.m(1,:),2);
                    po.m   = po.m(:,ix);
                    pr.m   = pr.m(:,ix);
                    po.b   = po.b(1,ix);
                    pr.b   = pr.b(1,ix);
                    po.n   = po.n(1,ix);
                    pr.n   = pr.n(1,ix);
                    po.W   = po.W(:,:,ix);
                    pr.W   = pr.W(:,:,ix);
                end
            else       
                % Use moments to compute means and variances, and then use these
                % to initialise the Gaussians
                %--------------------------------------------------

                [lq,~,mn,~,vr] = spm_GaussiansFromSuffStats(vb,mom,vr0);
                
                if logtpm.uniform
                    % Sort according to m values of first channel
                    [~,ix] = sort(mn(1,:),2);
                    mn     = mn(:,ix);
                end
            end 
            
            if isfield(obj,'mg') && (isfield(obj,'po') && isfield(obj,'pr') || isfield(obj,'mn') && isfield(obj,'vr'))
                % Parameters known
                %----------------------------------------------------------
                mg = obj.mg;
                if vb
                    po = obj.po;
                    pr = obj.pr;
                    
                    lq = spm_GaussiansFromSuffStats(vb,mom,po,pr);  
                else
                    mn = obj.mn;
                    vr = obj.vr;
                    
                    lq = spm_GaussiansFromSuffStats(vb,mom,vr0,mn,vr);
                end
            end
        else
            % Starting estimates for histograms
            %-----------------------------------------------------------------------
            lq = 0;
            
            for n=1:N
                maxval = -Inf;
                minval =  Inf;
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    maxval = max(max(buf(z).f{n}),maxval);
                    minval = min(min(buf(z).f{n}),minval);
                end
                maxval = max(maxval*1.5,-minval*0.05); % Account for bias correction effects
                minval = min(minval*1.5,-maxval*0.05);
                chan(n).interscal = [1 minval; 1 maxval]\[1;K];
                h0     = zeros(K,Kb);
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    cr       = round(buf(z).f{n}.*buf(z).bf{n}*chan(n).interscal(2) + chan(n).interscal(1));
                    cr       = min(max(cr,1),K);
                    for k1=1:Kb
                        h0(:,k1) = h0(:,k1) + accumarray(cr,buf(z).mu(:,k1),[K,1]);
                    end
                end
                chan(n).hist = h0;
            end
        end
    end

    for iter1=1:niter1
        if use_mog
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate cluster parameters
            %------------------------------------------------------------            
           
            for subitmog=1:nsubitmog                      
                oll = ll;                                              
                ll  = llr + llrb + lq; 
                mgm = 0;    
                mom = mom_struct(K,N,nomiss); 
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    
                    % Get template
                    mu  = double(buf(z).mu);
                    
                    % Denominator part of weight (wp) update
                    s   = 1./(mu*wp');                                        
                    mgm = mgm + s'*mu; 
                                        
                    % Evaluate the responsibilities using the current
                    % parameters values
                    [q,dll] = latent(buf(z).f,buf(z).bf,mg,mn,vr,po,mu,lkp,wp,buf(z).msk,buf(z).code,vb);                    
                    ll      = ll + dll;                                        
                
                    % Calculate sufficient statistics
                    mom = suffstats(mom,buf(z),q);
                end                            
                clear mu s q                                 
                
                % Plot/print objective value
                my_fprintf('MOG:\t%g\t%g\t%g\t%g\n',ll,lq,llr,llrb,verbose);
                if subitmog>1 || iter>1
                    spm_plot_convergence('Set',ll);
                end
                          
                % Restimate parameters using the current responsibilities
                if vb
                    % Posteriors from moments
                    [lq,s0,po] = spm_GaussiansFromSuffStats(vb,mom,po,pr);                    
                else
                    % Means and Variances from moments
                    [lq,s0,mn,vr] = spm_GaussiansFromSuffStats(vb,mom,vr0,mn,vr);
                end

                % Update within tissue mixing proportions
                for k=1:K
                    tmp   = s0(lkp==lkp(k));
                    mg(k) = s0(k)/sum(tmp);  % US eq. 27 (partly)
                end
                
                % Update tissue weights
                for k1=1:Kb
                    wp(k1) = (sum(s0(lkp==k1)) + wp_reg*1)/(mgm(k1) + wp_reg*Kb); % bias the solution towards 1
                end
                wp = wp/sum(wp); 

                if subitmog>1 && abs(ll-oll)<tol1*nm
                    % Improvement is small, so go to next step
                    break;
                end
            end  
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate histogram parameters
            %------------------------------------------------------------

            % Compute regularisation for histogram smoothing
            for n=1:N
               %x = (1:K)';
                for k1=1:Kb
                   %mom0 = sum(chan(n).hist(:,k1)) + eps;
                   %mom1 = sum(chan(n).hist(:,k1).*x) + eps;
                   %chan(n).lam(k1) = sum(chan(n).hist(:,k1).*(x-mom1./mom0).^2+1)/(mom0+1)+1;
                    chan(n).lam(k1) = Kb^2*double(vr0(N,N)*chan(n).interscal(2)^2);
                end
            end

            for subith=1:20
                oll  = ll;
                ll   = llr + llrb;
                for n=1:N
                    chan(n).lik  = spm_smohist(chan(n).hist,chan(n).lam);
                    chan(n).lik  = chan(n).lik*chan(n).interscal(2);
                    chan(n).alph = log(chan(n).lik+eps);
                    chan(n).hist = zeros(K,Kb);
                end
                mgm  = zeros(1,Kb);
                for z=1:length(z0)
                    mu   = double(buf(z).mu);
                    s   = 1./(mu*wp');
                    mgm = mgm + s'*mu;

                    [q,dll] = latent_nonpar(buf(z).f,buf(z).bf,chan,buf(z).mu,wp);
                    ll      = ll + dll;

                    cr  = cell(N,1);
                    for n=1:N
                        tmp   = buf(z).f{n}.*buf(z).bf{n}*chan(n).interscal(2) + chan(n).interscal(1);
                        cr{n} = min(max(round(tmp),1),K);
                    end
                    clear tmp
                    for k1=1:Kb
                        for n=1:N
                            chan(n).hist(:,k1) = chan(n).hist(:,k1) + accumarray(cr{n},q(:,k1),[K,1]);
                        end
                    end
                    clear cr
                end
                wp = (sum(chan(1).hist)+wp_reg*1)./(mgm+wp_reg*Kb);
                wp = wp/sum(wp);

                my_fprintf('Hist:\t%g\t%g\t%g\n',ll,llr,llrb,verbose);

                if subith>1 || iter>1
                    spm_plot_convergence('Set',ll); 
                end
                if subith>1 && abs(ll-oll)<tol1*nm
                    % Improvement is small, so go to next step
                    break;
                end
            end
            for n=1:N
                chan(n).lik  = spm_smohist(chan(n).hist,chan(n).lam);
                chan(n).lik  = chan(n).lik*chan(n).interscal(2);
                chan(n).alph = log(chan(n).lik+eps);
                chan(n).grad1 = convn(chan(n).alph,[0.5 0 -0.5]'*chan(n).interscal(2),  'same');
                chan(n).grad2 = convn(chan(n).alph,[1  -2  1  ]'*chan(n).interscal(2)^2,'same');
            end
        end
                
        if iter1 > 1 && ~(abs(ll-ooll)>2*tol1*nm)
            break; 
        end
        ooll = ll;

        if dobias
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate bias
            % Note that for multi-spectral data, the covariances among
            % channels are not computed as part of the second derivatives.
            % The aim is to save memory, and maybe make the computations
            % faster.
            %------------------------------------------------------------

            % Get means and precisions
            pr_bf = zeros(N,N,K); 
            if use_mog && vb
                for k=1:K, 
                    pr_bf(:,:,k) = po.n(k)*po.W(:,:,k); 
                end
                mn_bf = po.m;
            elseif use_mog && ~vb
                for k=1:K, 
                    pr_bf(:,:,k) = inv(vr(:,:,k)); 
                end
                mn_bf = mn;
            end
            
            for subitbf=1:nsubitbf
                for n=1:N
                    d3  = numel(chan(n).T);
                    if d3>0
                        % Compute objective function and its 1st and second derivatives
                        Alpha = zeros(d3,d3); % Second derivatives
                        Beta  = zeros(d3,1);  % First derivatives
                        ll    = llr + llrb + lq;
                        for z=1:length(z0)
                            if ~buf(z).nm, continue; end

                            if use_mog 
                                [q,dll] = latent(buf(z).f,buf(z).bf,mg,mn,vr,po,double(buf(z).mu),lkp,wp,buf(z).msk,buf(z).code,vb);
                                ll      = ll + dll;
                                
                                cr = cell(N,1);
                                for n1=1:N, cr{n1} = double(buf(z).f{n1}).*double(buf(z).bf{n1}); end

                                w1 = zeros(buf(z).nm,1);
                                w2 = zeros(buf(z).nm,1);
                                for k=1:K
                                    qk  = q(:,k);
                                    w0  = zeros(buf(z).nm,1);
                                    for n1=1:N
                                        w0 = w0 + pr_bf(n1,n,k)*(mn_bf(n1,k) - cr{n1});
                                    end
                                    w1  = w1 + qk.*w0;
                                    w2  = w2 + qk*pr_bf(n,n,k);
                                end
                                wt1   = zeros(d(1:2));
                                wt1(buf(z).msk{1}) = -(1 + cr{n}.*w1); % US eq. 34 (gradient)
                                wt2   = zeros(d(1:2));
                                wt2(buf(z).msk{1}) = cr{n}.*cr{n}.*w2 + 1; % Simplified Hessian of US eq. 34
                                clear cr
                            else
                                [q,dll] = latent_nonpar(buf(z).f,buf(z).bf,chan,buf(z).dat,wp);
                                ll      = ll + dll;
                                
                                cr0 = buf(z).f{n}.*buf(z).bf{n};
                                cr  = cr0*chan(n).interscal(2) + chan(n).interscal(1);
                                cr  = min(max(round(cr),1),K);
                                wt1 = zeros(d(1:2)); 
                                wt2 = zeros(d(1:2));
                                for k1=1:Kb
                                    qk  = q(:,k1);
                                    gr1 = chan(n).grad1(:,k1);
                                    gr1 = gr1(cr);
                                    gr2 = chan(n).grad2(:,k1);
                                    gr2 = min(gr2(cr),0); % Regularise
                                    wt1(buf(z).msk{1}) = wt1(buf(z).msk{1}) - qk.*(gr1.*cr0 + 1);
                                   %wt2(buf(z).msk{1}) = wt2(buf(z).msk{1}) - qk.*(gr1.*cr0 + gr2.*cr0.^2);
                                    wt2(buf(z).msk{1}) = wt2(buf(z).msk{1}) + qk.*(1 - gr2.*cr0.^2);
                                end
                            end

                            b3    = chan(n).B3(z,:)';
                            Beta  = Beta  + kron(b3,spm_krutil(wt1,chan(n).B1,chan(n).B2,0));
                            Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,chan(n).B1,chan(n).B2,1));
                            clear wt1 wt2 b3
                        end

                        oll     = ll;
                        C       = chan(n).C; % Inverse covariance of priors
                        oldT    = chan(n).T;

                        % Gauss-Newton update of bias field parameters
                        Update  = reshape((Alpha + C)\(Beta + C*chan(n).T(:)),size(chan(n).T));
                        clear Alpha Beta

                        % Line search to ensure objective function improves
                        armijo = 1.0;
                        for line_search=1:12
                            chan(n).T = chan(n).T - armijo*Update; % Backtrack if necessary

                            % Re-generate bias field, and compute terms of the objective function
                            chan(n).ll = double(-0.5*chan(n).T(:)'*C*chan(n).T(:));
                            for z=1:length(z0)
                                if ~buf(z).nm, continue; end
                                
                                bf             = transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T);
                                tmp            = bf(buf(z).msk{1});
                                if ~vb
                                    chan(n).ll = chan(n).ll + double(sum(tmp));
                                end
                                buf(z).bf{n}   = single(exp(tmp));
                            end
                            
                            llrb = 0;
                            for n1=1:N, 
                                llrb = llrb + chan(n1).ll; 
                            end
                            ll   = llr + llrb;
                            
                            mom = mom_struct(K,N,nomiss); 
                            for z=1:length(z0)
                                if ~buf(z).nm, continue; end
                                
                                if use_mog
                                    % Calculate responsibilities
                                    [q,dll] = latent(buf(z).f,buf(z).bf,mg,mn,vr,po,double(buf(z).mu),lkp,wp,buf(z).msk,buf(z).code,vb);
                    
                                    % Calculate moments
                                    mom = suffstats(mom,buf(z),q);
                                    clear q  
                                else
                                    [~,dll] = latent_nonpar(buf(z).f,buf(z).bf,chan,double(buf(z).mu),wp);                                                                    
                                end                                
                                ll = ll + dll;
                            end 
                            
                            if use_mog
                               olq = lq;
                               if vb                               
                                    lq = spm_GaussiansFromSuffStats(vb,mom,po,pr);                    
                               else                                    
                                    lq = spm_GaussiansFromSuffStats(vb,mom,vr0,mn,vr);
                               end    
                               ll = ll + lq;
                            end
                            
                            if ll>=oll
                                spm_plot_convergence('Set',ll); 
                                my_fprintf('Bias-%d:\t%g\t%g\t%g\t%g :o)\n',n,ll,lq,llr,llrb,verbose);
                                break;
                            else
                                my_fprintf('Bias-%d:\t%g\t%g\t%g\t%g :o(\n',n,ll,lq,llr,llrb,verbose);
                                ll        = oll;
                                lq        = olq;
                                chan(n).T = oldT;
                                armijo    = armijo*0.5;                                
                            end 
                        end
                        clear oldT
                    end
                end
                if subitbf > 1 && ~(ll-oll>tol1*nm)
                    % Improvement is only small, so go to next step
                    break;
                end
            end            
        end       

        if ~isempty(fig{1}) && use_mog   
            % Visualise responsibilities  
            q  = latent(buf(zix).f,buf(zix).bf,mg,mn,vr,po,double(buf(zix).mu),lkp,wp,buf(zix).msk,buf(zix).code,vb);            
            K1 = floor(sqrt(K));
            K2 = ceil(K/K1); 
            set(0,'CurrentFigure',fig{1});                    
            for i=1:K      
                subplot(K1,K2,i);
                slice(buf(zix).msk{1}) = q(:,i);
                imagesc(slice'); axis image xy off; title(['k=' num2str(lkp(i))]); colormap(gray);
            end 
            clear q
            drawnow
        end

        if ~isempty(fig{3})                
            % Visualise bias field
            set(0,'CurrentFigure',fig{3});                    
            for i=1:N
                subplot(1,N,i);
                slice(buf(zix).msk{i}) = buf(zix).bf{i};
                imagesc(slice'); axis image xy off; colormap(gray);
            end  
            drawnow
        end
            
        if iter==1 && iter1==1
            % Most of the log-likelihood improvements are in the first iteration.
            % Show only improvements after this, as they are more clearly visible.
            spm_plot_convergence('Clear'); 
            spm_plot_convergence('Init','Processing','Log-likelihood','Iteration'); 

            if numel(obj.lkp) ~= numel(lkp)
                if vb
                    [Kb,K,lkp,mg,po,pr] = mog_in_mog(vb,obj.lkp,po,pr);
                else
                    [Kb,K,lkp,mg,mn,vr] = mog_in_mog(vb,obj.lkp,mn,vr);
                end
            end
        end
    end

    if dodef
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate deformations
        %------------------------------------------------------------        
        
        % Compute likelihoods, and save them in buf.mu
        ll_const = 0;
        ll       = llr + llrb + lq;  
        if use_mog
            if vb
                % VB-MoG
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end                  
                    [~,dll,qt] = latent(buf(z).f,buf(z).bf,mg,mn,vr,po,double(buf(z).mu),lkp,wp,buf(z).msk,buf(z).code,vb);
                    ll         = ll + dll;
                    max_qt     = max(qt,[],2);                    
                    q = zeros([buf(z).nm Kb]);
                    for k1=1:Kb
                        for k=find(lkp==k1)
                            msk = qt(:,k)~=0;
                            q(msk,k1) = q(msk,k1) + exp(qt(msk,k) - max_qt(msk));
                        end
                        buf(z).mu(:,k1) = single(q(:,k1));
                    end                                    
                end            
            else
                % Regular MoG
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    qt       = log_likelihoods(buf(z).f,buf(z).bf,mg,mn,vr,po,buf(z).msk,buf(z).code,vb);                     
                    max_qt   = max(qt,[],2);
                    ll_const = ll_const + sum(max_qt);  
                    mu       = double(buf(z).mu);
                    mu       = bsxfun(@times,mu,wp);
                    mu       = bsxfun(@times,mu,1./sum(mu,2));
                    q = zeros([buf(z).nm Kb]);
                    for k1=1:Kb
                        for k=find(lkp==k1)
                            msk = qt(:,k)~=0;
                            q(msk,k1) = q(msk,k1) + exp(qt(msk,k) - max_qt(msk));
                        end
                        buf(z).mu(:,k1) = single(q(:,k1));
                    end                                

                    ll = ll + sum(log(sum(q.*mu + tiny,2)));
                end                
                clear mu
            end
            clear qt max_qt q
        else
            for z=1:length(z0)
                if ~buf(z).nm, continue; end
                q        = log_likelihoods_nonpar(buf(z).f,buf(z).bf,chan);
                max_q    = max(q,[],2);
                ll_const = ll_const + sum(max_q);
                q        = exp(bsxfun(@minus,q,max_q));
                mu       = bsxfun(@times,double(buf(z).mu),wp);
                mu       = bsxfun(@times,mu,1./sum(mu,2));
                ll       = ll + sum(log(sum(q.*mu + tiny,2)),1);
                buf(z).mu = single(q);
            end            
        end
        ll = ll + ll_const;
        
        oll = ll;        
        for subitdef=1:nitdef
            Alpha  = zeros([size(x0),numel(z0),6],'single');
            Beta   = zeros([size(x0),numel(z0),3],'single');
            for z=1:length(z0)
                if ~buf(z).nm, continue; end

                % Deformations from parameters
                [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,buf(z).msk);

                % Tissue probability map and spatial derivatives
                [mu,db1,db2,db3] = spm_sample_logpriors8(logtpm,x1,y1,z1);
                clear x1 y1 z1

                % Adjust for tissue weights
                s   = zeros(size(mu{1}));
                ds1 = zeros(size(mu{1}));
                ds2 = zeros(size(mu{1}));
                ds3 = zeros(size(mu{1}));
                for k1=1:Kb
                    mu{k1}   = wp(k1)*mu{k1};
                    db1{k1} = wp(k1)*db1{k1};
                    db2{k1} = wp(k1)*db2{k1};
                    db3{k1} = wp(k1)*db3{k1};
                    s       =  s  + mu{k1};
                    ds1     = ds1 + db1{k1};
                    ds2     = ds2 + db2{k1};
                    ds3     = ds3 + db3{k1};
                end
                for k1=1:Kb
                    mu{k1}   = mu{k1}./s;
                    db1{k1} = (db1{k1}-mu{k1}.*ds1)./s;
                    db2{k1} = (db2{k1}-mu{k1}.*ds2)./s;
                    db3{k1} = (db3{k1}-mu{k1}.*ds3)./s;
                end
                clear s ds1 ds2 ds3

                % Rotate gradients (according to initial affine registration) and
                % compute the sums of the tpm and its gradients, times the likelihoods
                % (from buf.mu).
                p   = zeros(buf(z).nm,1) + eps;
                dp1 = zeros(buf(z).nm,1);
                dp2 = zeros(buf(z).nm,1);
                dp3 = zeros(buf(z).nm,1);
                MM  = M*MT; % Map from sampled voxels to atlas data
                for k1=1:Kb
                    pp  = double(buf(z).mu(:,k1));
                    p   = p   + pp.*mu{k1};
                    dp1 = dp1 + pp.*(MM(1,1)*db1{k1} + MM(2,1)*db2{k1} + MM(3,1)*db3{k1});
                    dp2 = dp2 + pp.*(MM(1,2)*db1{k1} + MM(2,2)*db2{k1} + MM(3,2)*db3{k1});
                    dp3 = dp3 + pp.*(MM(1,3)*db1{k1} + MM(2,3)*db2{k1} + MM(3,3)*db3{k1});
                end
                clear mu db1 db2 db3

                % Compute first and second derivatives of the matching term.  Note that
                % these can be represented by a vector and tensor field respectively.
                tmp             = zeros(d(1:2));
                tmp(buf(z).msk{1}) = dp1./p; dp1 = tmp;
                tmp(buf(z).msk{1}) = dp2./p; dp2 = tmp;
                tmp(buf(z).msk{1}) = dp3./p; dp3 = tmp;
                
                Beta(:,:,z,1)   = -dp1;     % First derivatives
                Beta(:,:,z,2)   = -dp2;
                Beta(:,:,z,3)   = -dp3;

                Alpha(:,:,z,1)  = dp1.*dp1; % Second derivatives
                Alpha(:,:,z,2)  = dp2.*dp2;
                Alpha(:,:,z,3)  = dp3.*dp3;
                Alpha(:,:,z,4)  = dp1.*dp2;
                Alpha(:,:,z,5)  = dp1.*dp3;
                Alpha(:,:,z,6)  = dp2.*dp3;
                clear tmp p dp1 dp2 dp3
            end

            % Heavy-to-light regularisation
            if ~isfield(obj,'Twarp')
                switch iter
                case 1
                    prm = [param(1:3) 256*param(4:8)];
                case 2
                    prm = [param(1:3) 128*param(4:8)];
                case 3
                    prm = [param(1:3)  64*param(4:8)];
                case 4
                    prm = [param(1:3)  32*param(4:8)];
                case 5
                    prm = [param(1:3)  16*param(4:8)];
                case 6
                    prm = [param(1:3)  8*param(4:8)];
                case 7
                    prm = [param(1:3)  4*param(4:8)];
                case 8
                    prm = [param(1:3)  2*param(4:8)];
                otherwise
                    prm = [param(1:3)    param(4:8)];
                end
            else
                prm = [param(1:3)   param(4:8)];
            end

            % Add in the first derivatives of the prior term
            Beta   = Beta  + spm_diffeo('vel2mom',bsxfun(@times,Twarp,1./sk4),prm);

            % Gauss-Newton increment
            Update = bsxfun(@times,spm_diffeo('fmg',Alpha,Beta,[prm 2 2]),sk4);

            % Line search to ensure objective function improves
            armijo = 1.0;
            for line_search=1:12
                Twarp1 = Twarp - armijo*Update; % Backtrack if necessary

                % Recompute objective function
                llr1 = -0.5*sum(sum(sum(sum(Twarp1.*bsxfun(@times,spm_diffeo('vel2mom',bsxfun(@times,Twarp1,1./sk4),prm),1./sk4)))));
                ll1  = llr1 + llrb;
                ll1  = ll1  + ll_const;
                mom  = mom_struct(K,N,nomiss);
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    
                    [x1,y1,z1] = defs(Twarp1,z,x0,y0,z0,M,buf(z).msk);
                    mu         = spm_sample_logpriors8(logtpm,x1,y1,z1);
                    clear x1 y1 z1
                    
                    if vb
                        nmu = zeros([numel(mu{1}) K]);
                        for k1=1:Kb
                           nmu(:,k1) = double(mu{k1});
                        end
                        mu = nmu;
                        clear nmu
                        
                        [q,dll] = latent(buf(z).f,buf(z).bf,mg,mn,vr,po,mu,lkp,wp,buf(z).msk,buf(z).code,vb);
                        ll1     = ll1 + dll;
                        
                        mom   = suffstats(mom,buf(z),q);    
                        clear q
                    else                    
                        s   = zeros(size(mu{1}));
                        for k1=1:Kb, mu{k1} = mu{k1}*wp(k1); s = s + mu{k1}; end
                        for k1=1:Kb, mu{k1} = mu{k1}./s; end

                        sq = zeros(buf(z).nm,1);
                        for k1=1:Kb
                            sq = sq + double(buf(z).mu(:,k1)).*double(mu{k1});
                        end
                        
                        ll1 = ll1 + sum(log(sq));
                        clear sq
                    end
                    clear mu
                end

                lq1 = lq;
                if vb                       
                    lq1  = spm_GaussiansFromSuffStats(vb,mom,po,pr);
                    ll1 = ll1 + lq1;
                end
            
                if ll1<ll
                    % Still not better, so keep searching inwards.                    
                    my_fprintf('Warp:\t%g\t%g\t%g\t%g :o(\t(%g)\n',ll1,lq,llr1,llrb,armijo,verbose);                    
                    armijo = armijo*0.75;
                else
                    % Better. Accept the new solution.
                    spm_plot_convergence('Set',ll1); 
                    my_fprintf('Warp:\t%g\t%g\t%g\t%g :o)\t(%g)\n',ll1,lq,llr1,llrb,armijo,verbose);
                    ll     = ll1;
                    lq     = lq1;
                    llr    = llr1;
                    Twarp  = Twarp1;
                    break
                end
            end
            clear Alpha Beta
            
            if ~(abs(ll-oll)>tol1*nm)
                % Registration no longer helping, so move on
                break
            end
            oll = ll;
        end
        
        if ~isempty(fig{2})                
            % Visualise deformed template
            [x1,y1,z1] = defs(Twarp,zix,x0,y0,z0,M,buf(zix).msk);
            mu         = spm_sample_logpriors8(logtpm,x1,y1,z1);
            clear x1 y1 z1

            K1 = floor(sqrt(Kb));
            K2 = ceil(Kb/K1); 
            set(0,'CurrentFigure',fig{2});                                        
            for i=1:Kb     
                subplot(K1,K2,i);
                slice(buf(zix).msk{1}) = mu{i};
                imagesc(slice'); axis image xy off; colormap(gray);
            end 
            clear mu
            drawnow
        end
    end

    if ~(abs(ll-ooll)>2*tol1*nm) && iter>9
        % Finished
        break
    end
end

if nargout>=2 && dotpm && use_mog         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For estimating a template using maximum likelihood
    %------------------------------------------------------------
    
    % Compute likelihoods
    %-----------------------------------------------------------
    pf = zeros([prod(d(1:2)) d(3) Kb],'single'); % likelihood (p(data|parameters))
    for z=1:length(z0)
        if ~buf(z).nm, continue; end
              
        qt     = log_likelihoods(buf(z).f,buf(z).bf,mg,mn,vr,po,buf(z).msk,buf(z).code,vb);
        q      = zeros([buf(z).nm Kb]);
        max_qt = max(qt,[],2);
        for k1=1:Kb
            for k=find(lkp==k1)
                q(:,k1) = q(:,k1) + exp(qt(:,k) - max_qt);
            end
            pf(buf(z).msk{n},z,k1) = single(q(:,k1));
        end
        clear qt max_qt q
    end          
    pf = reshape(pf,[d Kb]);    
    clear buf

    if ~isempty(fig{4}) 
        % Display likelihoods in subject space
        set(0,'CurrentFigure',fig{4});      
        for i=1:Kb
            subplot(2,Kb,i);
            imagesc(pf(:,:,zix,i)'); axis image xy off; colormap(gray);
        end  
    end
    
    % Push likelihoods (in subject space) to template space 
    %-----------------------------------------------------------
    t       = identity(d) + bsxfun(@times,Twarp,1./sk4); 
    MM      = M*MT;
    t       = affine_transf(MM,t);         
    [pf,c] = spm_diffeo('pushc',pf,t,dtpm);
    clear t                
    
    msk       = isfinite(pf);    
    msk       = bsxfun(@and,msk,c>0);
    pf(~msk) = 0;
    clear msk c
    
    pf = max(pf,eps);
    pf = bsxfun(@rdivide,pf,sum(pf,4));    
    
    if ~isempty(fig{4}) 
        % Display likelihoods pushed to template space
        set(0,'CurrentFigure',fig{4});      
        for i=(Kb + 1):2*Kb
            subplot(2,Kb,i);
            imagesc(pf(:,:,floor(logtpm.d(3)/2 + 1),i - Kb)'); axis image xy off; colormap(gray);
        end  
    end
    
    % Solve dL/dmu=0 for mu (where L is the objective function and mu is the template)  
    %-----------------------------------------------------------
    
    nmu = prod(dtpm); % number of template voxels
    
    % log(mu)
    logmu = zeros([Kb nmu],'single'); 
    for k=1:Kb
        logmu(k,:) = reshape(logtpm.dat{k},[1 nmu]);
    end
    clear logtpm
        
    % exp(log(w) + log(mu)) = w*mu
    logwp = repmat(log(wp)',1,nmu);    
    wmu   = bsxfun(@plus,logmu,logwp);    
    wmu   = exp(wmu);
    clear logmu

    % w*mu*P(f|theta)    
    pf    = reshape(pf,[nmu Kb])'; 
    wmupf = bsxfun(@times,wmu,pf);
    clear pf

    % w*mu*P(f|theta)/sum(w*mu*P(f|theta),k)
    munum = bsxfun(@rdivide,wmupf,sum(double(wmupf),1)); % numerator of update
    munum(~isfinite(munum)) = 0;
    clear wmupf
    
    % w/sum(w*mu,k)
    wmu   = 1./sum(double(wmu),1);     
    muden = bsxfun(@times,wmu,exp(logwp)); % denominator of update
    muden(~isfinite(muden)) = 0;    
else
    munum = 0;
    muden = 0;
end

% Save the results
obj.Twarp = Twarp;
obj.Tbias = {chan(:).T};
obj.mg    = mg;
obj.wp    = wp;
obj.ll    = ll;
obj.nm    = nm;
if use_mog && vb
    obj.po = po;
    obj.pr = pr;
elseif use_mog
    obj.mn = mn;
    obj.vr = vr;
end
    
return;
%=======================================================================

%=======================================================================
function [Q,ll,L] = latent(f,bf,mg,mn,vr,po,mu,lkp,wp,msk,code,vb)
logmu = log_spatial_priors(mu,wp);
L     = log_likelihoods(f,bf,mg,mn,vr,po,msk,code,vb);

N  = numel(f);
K  = numel(lkp);
Kb = max(lkp);
% for i=2:2^N 
%     msk0 = dec2bin(i-1,N)=='1';  
%     ind  = find(code==msk0*(2.^(0:(N-1))'));
%     if ~isempty(ind)
%         for k1=1:Kb
%             for k=find(lkp==k1)
%                 Q(ind,k) = Q(ind,k) + logmu(ind,k1);
%             end
%         end
%     end
% end  
Q = zeros(size(L));
for k1=1:Kb
    for k=find(lkp==k1)
        Q(:,k) = L(:,k) + logmu(:,k1);
    end
end

if vb
    logSumQ = logsumexp(Q,2);
    logQ    = bsxfun(@minus,Q,logSumQ);
    Q       = exp(logQ);
    clear logSumQ
    
    % Compute lower bound
    %----------------------------------------------------------------------
    
    % Responsibilities
    ll = -sum(sum(Q.*logQ));
    clear logQ

    % Bias field    
    nbf = ones([nnz(msk{1}) N]);
    for n=1:N
        nbf(:,n) = double(bf{n});
    end
    bf  = nbf;
    clear nbf
    bf  = repmat(log(prod(bf,2)),1,K); % Log|bf|
    ll  = ll + sum(sum(Q.*bf));
    clear bf
    
    % Template
    if size(Q,2)~=size(mu,2)        
%         mu  = bsxfun(@times,mu(:,lkp),mg);        
%         for k=1:Kb
%             kk  = find(lkp==k);
%             smu = sum(mu(:,kk),2);
%             for k1=kk                
%                 mu(:,k1) = mu(:,k1)./smu;            
%             end
%         end
%         mu    = bsxfun(@times,mu,wp(lkp)); 
%         logmu = log(mu);
        mu    = bsxfun(@times,mu,wp);
        mu    = bsxfun(@times,mu,1./sum(mu,2));
        mu    = bsxfun(@times,mu(:,lkp),mg);
        logmu = log(mu);
    end
    ll = ll + sum(sum(Q.*logmu));
else
    [Q,ll] = safe_softmax(Q);
end
%=======================================================================

%=======================================================================
function L = log_likelihoods(f,bf,mg,mn,vr,po,msk,code,vb)
K = numel(mg);
N = numel(f);
d = size(msk{1});

cr = zeros([prod(d) N]);
for n=1:N
    cr(msk{n},n) = double(f{n}).*double(bf{n});
end

if vb
    nbf = ones([prod(d) N]);
    for n=1:N
        nbf(msk{n},n) = double(bf{n});
    end
    bf  = nbf;
    clear nbf
end

L = zeros([prod(d) K]);
for i=2:2^N    
    msk0 = dec2bin(i - 1,N)=='1';  
    ind  = find(code==msk0*(2.^(0:(N-1))'));
    
    if ~isempty(ind)
        for k=1:K
            if vb  
                t1         = psi(0, 0.5*repmat(po.n(k)+1,N,1) - 0.5*[1:N]');
                lnLamTilde = sum(t1) + N*log(2)  + LogDet(po.W(msk0,msk0,k));    

                diff1      = bsxfun(@minus,cr(ind,msk0)',po.m(msk0,k));
                Q          = chol(po.W(msk0,msk0,k))*diff1;
                E          = N/po.b(k) + po.n(k)*dot(Q,Q,1);

                L(ind,k)   = 0.5*(lnLamTilde - E') + log(mg(k)) + log(prod(bf(ind,msk0),2)) - N/2*log(2*pi);
            else
                C          = chol(vr(msk0,msk0,k));
                diff1      = bsxfun(@minus,cr(ind,msk0),mn(msk0,k)')/C;
                L(ind,k)   = log(mg(k)) - (N/2)*log(2*pi) - sum(log(diag(C))) - 0.5*sum(diff1.*diff1,2);
            end
        end
    end
end

nL = zeros([nnz(L(:,1)) K]);
for k=1:K
   nL(:,k) = L(msk{1},k);
end
L = nL;
%=======================================================================

%=======================================================================
function mu = log_spatial_priors(mu,wp)
mu = bsxfun(@times,mu,wp);
mu = bsxfun(@times,mu,1./sum(mu,2));
mu = log(mu);
%=======================================================================

%=======================================================================
function L = log_likelihoods_nonpar(f,bf,chan)
K  = size(chan(1).lik,1);
Kb = size(chan(1).lik,2);
N  = numel(chan);
L  = zeros(numel(f{1}),Kb);
for n=1:N
    if isempty(bf)
        tmp = f{n}*chan(n).interscal(2) + chan(n).interscal(1);
    else
        tmp = f{n}.*bf{n}*chan(n).interscal(2) + chan(n).interscal(1);
    end
    tmp     = min(max(round(tmp),1),K);
    loglik  = chan(n).alph;
    for k1=1:Kb
        L(:,k1) = L(:,k1)+loglik(tmp,k1);
    end
end
%=======================================================================

%=======================================================================
function [Q,ll] = latent_nonpar(f,bf,chan,mu,wp)
mu     = log_spatial_priors(mu,wp);
Q      = log_likelihoods_nonpar(f,bf,chan);
Q      = Q + mu;
[Q,ll] = safe_softmax(Q);
%=======================================================================
      
%=======================================================================
function [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,msk)
x1a = x0    + double(Twarp(:,:,z,1));
y1a = y0    + double(Twarp(:,:,z,2));
z1a = z0(z) + double(Twarp(:,:,z,3));
if nargin>=7  
    x1a = x1a(msk{1});
    y1a = y1a(msk{1});
    z1a = z1a(msk{1});
end
[x1,y1,z1] = affine_transf(M,x1a,y1a,z1a);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%=======================================================================

%=======================================================================
function mom = suffstats(mom,buf,q)
d = size(buf.msk{1});
N = numel(buf.f);
K = numel(mom(1).s0);

X = zeros([prod(d(1:2)) N]);
for n=1:N
    X(buf.msk{n},n) = double(buf.f{n}).*double(buf.bf{n});    
end

Q = zeros([prod(d(1:2)) K]);
for k=1:K
    Q(buf.msk{n},k) = double(q(:,k));    
end

mom = spm_SuffStats(X,Q,buf.code,mom);  
%=======================================================================                    

%=======================================================================
function varargout = mog_in_mog(varargin)
% A crude heuristic to replace a single Gaussian by a bunch of Gaussians
% If there is only one Gaussian, then it should be the same as the
% original distribution.
vb  = varargin{1};
lkp = varargin{2};
K   = numel(lkp);
Kb  = max(lkp);
mg  = ones(1,K)/K;

if vb           
    % Posterior
    po1 = varargin{3};    
    
    m1 = po1.m;
    b1 = po1.b;
    W1 = po1.W;
    n1 = po1.n;   

    po.m = zeros(size(m1));
    po.b = zeros(size(b1));
    po.W = zeros(size(W1));
    po.n = zeros(size(n1));

    % Prior
    pr1 = varargin{4};
    
    m0 = pr1.m;
    b0 = pr1.b;
    W0 = pr1.W;
    n0 = pr1.n;
    
    pr.m = zeros(size(m0));
    pr.b = zeros(size(b0));
    pr.W = zeros(size(W0));
    pr.n = zeros(size(n0));
    
    N = size(m1,1);
else
    % Mean
    mn1 = varargin{3};
    mn  = ones(size(mn1));
    
    % Variance
    vr1 = varargin{4};        
    vr  = zeros(size(vr1));
    
    N = size(mn,1);
end

for k1=1:Kb
    kk 	        = sum(lkp==k1);
    mg(lkp==k1) = 1/kk;
    w  		    = 1./(1 + exp(-(kk-1)*0.25)) - 0.5; 

    if vb
        vr1               = inv(n1(k1)*W1(:,:,k1));
        pr1               = inv(vr1*(1 - w));                        
        po.b(lkp==k1)     = b1(k1);
        po.m(:,lkp==k1)   = sqrtm(vr1)*randn(N,kk)*w + repmat(m1(:,k1),[1,kk]);
        po.n(lkp==k1)     = n1(k1);
        po.W(:,:,lkp==k1) = repmat(pr1/n1(k1),[1 1 kk]);

        pr.b(lkp==k1)     = b0(k1);
        pr.m(:,lkp==k1)   = repmat(m0(:,k1),[1,kk]);
        pr.n(lkp==k1)     = n0(k1);
        pr.W(:,:,lkp==k1) = repmat(W0(:,:,k1),[1 1 kk]);        
    else
        mn(:,lkp==k1)   = sqrtm(vr1(:,:,k1))*randn(N,kk)*w + repmat(mn1(:,k1),[1,kk]);
        vr(:,:,lkp==k1) = repmat(vr1(:,:,k1)*(1 - w),[1,1,kk]);
    end                                        
end  

varargout{1} = Kb;
varargout{2} = K;
varargout{3} = lkp;
varargout{4} = mg;
if vb
    varargout{5} = po;
    varargout{6} = pr;
else
    varargout{5} = mn;
    varargout{6} = vr;
end
%=======================================================================

%=======================================================================
function count = my_fprintf(varargin)
if varargin{end}
    count = fprintf(varargin{1:end - 1});
else
    count = 0;
end
%=======================================================================