function obj = spm_preprocx(obj,logtpm,fig)
% New and improved (?) unified segmentation routine
%
% FORMAT obj = spm_preprocx(obj,logtpm)
%
%_______________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin<3, fig = cell(4,1); end

V         = obj.image;
descrip   = obj.descrip;
N         = numel(V);
Affine    = obj.Affine;
M         = logtpm.M\Affine*V(1).mat;
d0        = V(1).dim(1:3);
dtpm      = logtpm.d;
vx        = vxsize(V(1).mat);
sk        = max([1 1 1],round(obj.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
wp_reg    = obj.wp_reg;
tol1      = obj.tolseg;
tiny      = eps*eps;
lkp       = obj.lkp;
K         = numel(obj.lkp);
Kb        = max(obj.lkp);

kron = @(a,b) spm_krutil(a,b);

niter      = obj.niter;
niter1     = obj.niter1;
nsubitmog  = obj.nsubitmog;
nsubitbf   = obj.nsubitbf;
nitdef     = obj.nitdef;

dobias = obj.dobias;
dowp   = obj.dowp;
dotpm  = obj.dotpm;

missing_data = obj.missing_data;

def_done0 = obj.def_done;
def_done  = 0;

uniform = obj.iter==1 && ~obj.use_tpm; % The TPMs given are uniform distributions (i.e., run k-means)

dodef = obj.dodef && ~uniform;

bb = obj.bb;

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

if obj.iter==1
    % Allocate deformations        
    create_nii(obj.pth_def,zeros([d 3],'single'),eye(4),'float32','def');      
end

Nii   = nifti(obj.pth_def);
Twarp = single(Nii.dat(:,:,:,:)); 
clear Nii
llr   = -0.5*sum(sum(sum(sum(Twarp.*bsxfun(@times,spm_diffeo('vel2mom',bsxfun(@times,Twarp,1./sk4),param),1./sk4)))));

% Initialise bias correction
%-----------------------------------------------------------------------
cl   = cell(N,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});

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

% Load the data
%-----------------------------------------------------------------------

% Total number of voxels  
nm = 0;

% Data type of images
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
buf = struct('f',cl,'dat',cl,'bf',cl,'code',cl,'msk',cl);
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
        buf(z).msk{n} = buf(z).msk{n} & get_msk(fz{n},descrip);
        buf(z).nm{n}  = nnz(buf(z).msk{n});
    end    
    
    if ~missing_data
        msk = true;
        for n=1:N
            msk = msk & buf(z).msk{n};
        end

        for n=1:N
            buf(z).msk{n} = msk;
            buf(z).nm{n}  = nnz(buf(z).msk{n});
        end
    end
    clear msk
                
    if isfield(obj,'msk') && ~isempty(obj.msk)
        % Exclude any voxels to be masked out
        msk = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        for n=1:N
            buf(z).msk{n} = buf(z).msk{n} & msk;
        end
    end
       
    % Sum up intensities of all voxels (used for later normalising image
    % intensities)
    a = zeros(N,1); % Intensity of voxels (for each channel)
    for n=1:N
        a(n) = sum(fz{n}(buf(z).msk{n}));
    end
    
    % For dealing with missing data
    code = zeros([prod(d(1:2)) 1],typ);
    for n=1:N
        code = bitor(code,bitshift(feval(cast,buf(z).msk{n}(:)),(n-1)));
    end
    buf(z).code = code;
    
    buf(z).Nm = nnz(code~=0);
    nm        = nm + buf(z).Nm; 
    
    for n=1:N                                   
        % Normalise image intensities (not if CT)            
        if strcmp(descrip,'CT')
            a(n) = 1.0;
        else
            a(n) = 512/(a(n)/buf(z).nm{n});
        end
       
        % Eliminate unwanted voxels
        if scrand(n)
            % The variances used in a GMM may not be ideal because of the
            % discretisation of image intensities (Data is an integer
            % type), so to prevent aliasing in the histogram, small random
            % values are added.  It's not elegant, but the alternative
            % would be too slow for practical use.
            buf(z).f{n}  = single(a(n)*fz{n}(buf(z).msk{n})+rand(buf(z).nm{n},1)*scrand(n)-scrand(n)/2);
        else
            buf(z).f{n}  = single(a(n)*fz{n}(buf(z).msk{n}));
        end    

        % Create moments to be used for constructing a variance prior
        s0(n) = s0(n) + buf(z).nm{n};
        s1(n) = s1(n) + sum(buf(z).f{n});
        S2(n) = S2(n) + sum(buf(z).f{n}.^2);
    end

    % Create a buffer for tissue probability info
    buf(z).dat = zeros([buf(z).Nm Kb],'single');
end

% Set-up MoG struct
mog = obj.mog;
vb  = mog.vb;

% Construct a ``Wishart-style prior''
% Mostly here to prevent some of the numerical instabilities that arise
% from computing the variance from the moments in a single pass.            
mog.vr0 = diag(S2./s0 - (s1./s0).^2)/K^2;
clear s0 s1 S2

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
        bf             = transf(B1,B2,B3(z,:),T);
        tmp            = bf(buf(z).msk{n});
        if ~vb
            chan(n).ll = chan(n).ll + double(sum(tmp));
        end
        buf(z).bf{n}   = single(exp(tmp));
    end    
    llrb = llrb + chan(n).ll;    
    clear B1 B2 B3 T C tmp
end

% For debugging
%-----------------------------------------------------------------------
zix     = floor(d(3)/2 + 1);
verbose = obj.verbose;
slice   = zeros(d(1:2),'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run algorithm
%------------------------------------------------------------

spm_plot_convergence('Init','Initialising','Log-likelihood','Iteration');

ll = -Inf;
lb = 0;
lh = 0;
for iter=1:niter
    % Load the warped prior probability images into the buffer
    %------------------------------------------------------------
    for z=1:length(z0)
        if ~buf(z).Nm, continue; end
        
        [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,buf(z).code);        
        tpm        = spm_sample_logpriors8(logtpm,x1,y1,z1);
        for k1=1:Kb
            buf(z).dat(:,k1) = tpm{k1};
        end         
        clear x1 y1 z1
                
        if ~isempty(fig{2}) && z==zix && iter==1
            % Visualise initial template
            K1 = floor(sqrt(Kb));
            K2 = ceil(Kb/K1); 
            set(0,'CurrentFigure',fig{2});                                        
            for i=1:Kb            
                subplot(K1,K2,i);
                slice(buf(z).code>0) = tpm{i};
                imagesc(slice'); axis image xy off; colormap(gray);
                title(['TPM, k=' num2str(i)]);
            end 
            drawnow
        end
        clear tpm
    end
           
    if iter==1
        % Get initial moments
        %----------------------------------------------------------
        K   = Kb;
        lkp = 1:Kb;
                  
        if isfield(obj,'wp')
            wp = obj.wp;
        else
            wp = ones(1,Kb)/Kb;
        end

        if isfield(obj,'mg')
            mg = obj.mg;
        else
            mg = ones(1,Kb);
        end

        mom = mom_from_buf(buf,Kb,N,z0,uniform);
    end
    
    for iter1=1:niter1        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate cluster parameters
        %------------------------------------------------------------            
        for subitmog=1:nsubitmog 
            % M-step   
            %----------------------------------------------------------
            [~,~,mog,s0] = spm_GaussiansFromSuffStats(mom,mog);                  

            if dowp
                % Update tissue mixing weights
                mgm = 0;
                for z=1:length(z0)        
                    tpm = double(buf(z).dat);
                    s   = 1./(tpm*wp');                                        
                    mgm = mgm + s'*tpm; 
                end
                for k=1:Kb
                    wp(k) = (sum(s0(lkp==k)) + wp_reg*1)/(mgm(k) + wp_reg*Kb);
                end
                wp = wp/sum(wp);    
            end

            for k=1:K
                % Update within tissue mixing weights
                tmp     = s0(lkp==lkp(k));
                mg(1,k) = s0(k)/sum(tmp);
            end
            
            % E-step        
            %----------------------------------------------------------
            oll = ll;
            ll  = llr + llrb;
            mom = mom_struct(K,N); 
            for z=1:length(z0)     
                if ~buf(z).Nm, continue; end

                cr = zeros(prod(d(1:2)),N);
                for n=1:N
                    cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n});
                end

                % Calculate responsibilities
                [q,dlq] = latent(buf(z),double(buf(z).dat),mg,mog,wp,lkp,cr);
                ll      = ll + sum(dlq);

                % Calculate sufficient statistics                
                mom = spm_SuffStats(cr,q,buf(z).code,mom);
            end     

            if vb
                [lb,lh] = spm_GaussiansFromSuffStats(mom,mog);  
                ll      = ll + lb + lh;                 
            end

            % Plot/print objective value
            my_fprintf('MOG:\t%g\t%g\t%g\t%g\t%g\n',ll,lb,lh,llr,llrb,verbose);
            if subitmog>1 || iter>1
                spm_plot_convergence('Set',ll);
            end

            if subitmog>1 && abs(ll-oll)<tol1*nm
                % Improvement is small, so go to next step
                break;
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
            if vb
                for k=1:K 
                    pr_bf(:,:,k) = mog.po.n(k)*mog.po.W(:,:,k); 
                end
                mn_bf = mog.po.m;
            else
                for k=1:K
                    pr_bf(:,:,k) = inv(mog.vr(:,:,k)); 
                end
                mn_bf = mog.mn;
            end
            
            for subitbf=1:nsubitbf
                for n=1:N
                    d3  = numel(chan(n).T);
                    if d3>0
                        % Compute objective function and its 1st and second derivatives
                        Alpha = zeros(d3,d3); % Second derivatives
                        Beta  = zeros(d3,1);  % First derivatives                        
                        for z=1:length(z0)
                            if ~buf(z).Nm, continue; end                           

                            cr = zeros(prod(d(1:2)),N);
                            for n1=1:N
                                cr(buf(z).msk{n1},n1) = double(buf(z).f{n1}).*double(buf(z).bf{n1});
                            end

                            q = latent(buf(z),double(buf(z).dat),mg,mog,wp,lkp,cr);
                            
                            w1 = zeros(buf(z).Nm,1);
                            w2 = zeros(buf(z).Nm,1);
                            for k=1:K
                                qk  = q(buf(z).msk{n},k);
                                w0  = zeros(buf(z).Nm,1);
                                for n1=1:N
                                    w0 = w0 + pr_bf(n1,n,k)*(mn_bf(n1,k) - cr(buf(z).msk{n1},n1));
                                end
                                w1  = w1 + qk.*w0;
                                w2  = w2 + qk*pr_bf(n,n,k);
                            end
                            wt1                = zeros(d(1:2));
                            wt1(buf(z).code>0) = -(1 + cr(buf(z).msk{n},n).*w1); % US eq. 34 (gradient)
                            wt2                = zeros(d(1:2));
                            wt2(buf(z).code>0) = cr(buf(z).msk{n},n).*cr(buf(z).msk{n},n).*w2 + 1; % Simplified Hessian of US eq. 34
                            clear cr

                            b3    = chan(n).B3(z,:)';
                            Beta  = Beta  + kron(b3,spm_krutil(wt1,chan(n).B1,chan(n).B2,0));
                            Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,chan(n).B1,chan(n).B2,1));
                            clear wt1 wt2 b3
                        end

                        oll     = ll;
                        olb     = lb;
                        olh     = lh;
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
                                if ~buf(z).Nm, continue; end
                                
                                bf             = transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T);
                                tmp            = bf(buf(z).msk{n});
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
                            
                            mom1 = mom_struct(K,N); 
                            for z=1:length(z0)     
                                if ~buf(z).Nm, continue; end

                                cr = zeros(prod(d(1:2)),N);
                                for n1=1:N
                                    cr(buf(z).msk{n1},n1) = double(buf(z).f{n1}).*double(buf(z).bf{n1});
                                end

                                % Calculate responsibilities
                                [q,dlq] = latent(buf(z),double(buf(z).dat),mg,mog,wp,lkp,cr);
                                ll      = ll + sum(dlq);

                                if vb
                                    % Calculate sufficient statistics                
                                    mom1 = spm_SuffStats(cr,q,buf(z).code,mom1);
                                end 
                            end     

                            if vb
                                [lb,lh] = spm_GaussiansFromSuffStats(mom1,mog);  
                                ll      = ll + lb + lh;  
                            end
                            
                            if ll>=oll
                                spm_plot_convergence('Set',ll); 
                                my_fprintf('Bias-%d:\t%g\t%g\t%g\t%g\t%g :o)\n',n,ll,lb,lh,llr,llrb,verbose);
                                break;
                            else
                                my_fprintf('Bias-%d:\t%g\t%g\t%g\t%g\t%g :o(\n',n,ll,lb,lh,llr,llrb,verbose);
                                ll        = oll;
                                lb        = olb;
                                lh        = olh;
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

        if ~isempty(fig{1})  
            % Visualise responsibilities  
            q  = latent(buf(zix),double(buf(zix).dat),mg,mog,wp,lkp);
            K1 = floor(sqrt(K + 1));
            K2 = ceil((K + 1)/K1); 
            set(0,'CurrentFigure',fig{1}); clf(fig{1});            
            for i=1:K + 1                
                subplot(K1,K2,i);
                if i==K + 1
                    slice(buf(zix).msk{1}) = buf(zix).f{1};
                    imagesc(slice'); axis image xy off; title('f1'); colormap(gray);
                else
                    slice = reshape(q(:,i),size(buf(zix).msk{1}));
                    imagesc(slice'); axis image xy off; title(['q, k=' num2str(lkp(i))]); colormap(gray);
                end                
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
                title(['bf, n=' num2str(i)]);
            end  
            drawnow
        end
            
        if iter==1 && iter1==1
            % Most of the log-likelihood improvements are in the first iteration.
            % Show only improvements after this, as they are more clearly visible.
            spm_plot_convergence('Clear'); 
            spm_plot_convergence('Init','Processing','Log-likelihood','Iteration'); 

            if numel(obj.lkp) ~= numel(lkp)            
                [Kb,K,lkp,mg,mog,mom] = split_mog(obj.lkp,mog,mom);
            end
        end
    end

    if dodef && iter<niter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate deformations
        %------------------------------------------------------------                        
        
        % Compute responsibilities (VB) or likelihoods (MAP), and save them in buf.dat
        ll_const = 0;
        ll       = llr + llrb + lb + lh;          
        if vb
            for z=1:length(z0)
                if ~buf(z).Nm, continue; end  

                [qt,dlq] = latent(buf(z),double(buf(z).dat),mg,mog,wp,lkp); 
                ll       = ll + sum(dlq);

                q   = zeros([buf(z).Nm Kb]);
                msk = buf(z).code>0;
                for k1=1:Kb
                    for k=find(lkp==k1)
                        q(:,k1) = q(:,k1) + qt(msk,k);
                    end
                    buf(z).dat(:,k1) = single(q(:,k1));
                end     
            end            
        else
            for z=1:length(z0)
                if ~buf(z).Nm, continue; end

                tpm = double(buf(z).dat);
                tpm = bsxfun(@times,tpm,wp);
                tpm = bsxfun(@times,tpm,1./sum(tpm,2));

                qt       = log_likelihoods(K,buf(z),mg,mog);
                max_qt   = max(qt,[],2);
                ll_const = ll_const + nansum(max_qt);  
                q        = zeros([buf(z).Nm Kb]);
                msk      = buf(z).code>0;
                for k1=1:Kb
                    for k=find(lkp==k1)
                        q(:,k1) = q(:,k1) + exp(qt(msk,k) - max_qt(msk));
                    end
                    buf(z).dat(:,k1) = single(q(:,k1));
                end                                

                ll  = ll + sum(log(sum(q.*tpm + tiny,2)));
            end                
            clear tpm
        end
        clear qt max_qt q
        ll = ll + ll_const;
        
        oll = ll;        
        for subitdef=1:nitdef
            
            % Compute gradients (Beta) and Hessians (Alpha)
            Alpha  = zeros([size(x0),numel(z0),6],'single');
            Beta   = zeros([size(x0),numel(z0),3],'single');
            for z=1:length(z0)
                if ~buf(z).Nm, continue; end

                % Deformations from parameters
                [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,buf(z).code);

                % Tissue probability map and spatial derivatives
                [tpm,db1,db2,db3] = spm_sample_logpriors8(logtpm,x1,y1,z1);
                clear x1 y1 z1

                % Adjust for tissue weights
                s   = zeros(size(tpm{1}));
                ds1 = zeros(size(tpm{1}));
                ds2 = zeros(size(tpm{1}));
                ds3 = zeros(size(tpm{1}));
                for k1=1:Kb
                    tpm{k1}   = wp(k1)*tpm{k1};
                    db1{k1} = wp(k1)*db1{k1};
                    db2{k1} = wp(k1)*db2{k1};
                    db3{k1} = wp(k1)*db3{k1};
                    s       =  s  + tpm{k1};
                    ds1     = ds1 + db1{k1};
                    ds2     = ds2 + db2{k1};
                    ds3     = ds3 + db3{k1};
                end
                for k1=1:Kb
                    tpm{k1}   = tpm{k1}./s;
                    db1{k1} = (db1{k1}-tpm{k1}.*ds1)./s;
                    db2{k1} = (db2{k1}-tpm{k1}.*ds2)./s;
                    db3{k1} = (db3{k1}-tpm{k1}.*ds3)./s;
                end
                clear s ds1 ds2 ds3

                % Rotate gradients (according to initial affine registration) and
                % compute the sums of the tpm and its gradients, times the
                % likelihoods/responsibilities
                p   = zeros(buf(z).Nm,1) + eps;
                dp1 = zeros(buf(z).Nm,1);
                dp2 = zeros(buf(z).Nm,1);
                dp3 = zeros(buf(z).Nm,1);
                MM  = M*MT; % Map from sampled voxels to atlas data
                for k1=1:Kb
                    pp  = double(buf(z).dat(:,k1));
                    p   = p   + pp.*tpm{k1};
                    dp1 = dp1 + pp.*(MM(1,1)*db1{k1} + MM(2,1)*db2{k1} + MM(3,1)*db3{k1});
                    dp2 = dp2 + pp.*(MM(1,2)*db1{k1} + MM(2,2)*db2{k1} + MM(3,2)*db3{k1});
                    dp3 = dp3 + pp.*(MM(1,3)*db1{k1} + MM(2,3)*db2{k1} + MM(3,3)*db3{k1});
                end
                clear tpm db1 db2 db3

                % Compute first and second derivatives of the matching term.  Note that
                % these can be represented by a vector and tensor field respectively.
                tmp                = zeros(d(1:2));
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
            if dotpm
                switch (obj.iter - 2)
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
                if ~def_done0
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
                ll1  = llr1 + llrb + ll_const;                                
                olb  = lb;
                olh  = lh;
                mom1 = mom_struct(K,N); 
                for z=1:length(z0)
                    if ~buf(z).Nm, continue; end
                    
                    [x1,y1,z1] = defs(Twarp1,z,x0,y0,z0,M,buf(z).code);
                    tpm        = spm_sample_logpriors8(logtpm,x1,y1,z1);
                    clear x1 y1 z1
                    
                    if vb                                              
                        cr = zeros(prod(d(1:2)),N);
                        for n1=1:N
                            cr(buf(z).msk{n1},n1) = double(buf(z).f{n1}).*double(buf(z).bf{n1});
                        end

                        B = zeros(buf(z).Nm,Kb);
                        for k1=1:Kb
                            B(:,k1) = double(tpm{k1});
                        end
                        
                        % Calculate responsibilities
                        [q,dlq] = latent(buf(z),B,mg,mog,wp,lkp,cr);
                        ll1     = ll1 + sum(dlq);
                        
                        % Calculate sufficient statistics                
                        mom1 = spm_SuffStats(cr,q,buf(z).code,mom1);                                
                    else                    
                        s   = zeros(size(tpm{1}));
                        for k1=1:Kb, tpm{k1} = tpm{k1}*wp(k1); s = s + tpm{k1}; end
                        for k1=1:Kb, tpm{k1} = tpm{k1}./s; end

                        sq = zeros(buf(z).Nm,1);
                        for k1=1:Kb
                            sq = sq + double(buf(z).dat(:,k1)).*double(tpm{k1});
                        end
                        
                        ll1 = ll1 + sum(log(sq));
                        clear sq
                    end
                    clear tpm
                end

                if vb
                    [lb,lh] = spm_GaussiansFromSuffStats(mom1,mog);  
                    ll1     = ll1 + lb + lh;  
                end
            
                if ll1<ll
                    % Still not better, so keep searching inwards.                    
                    my_fprintf('Warp:\t%g\t%g\t%g\t%g\t%g :o(\t(%g)\n',ll1,lb,lh,llr1,llrb,armijo,verbose);                    
                    armijo = armijo*0.75;
                    lb     = olb;
                    lh     = olh;
                else
                    % Better. Accept the new solution.                    
                    spm_plot_convergence('Set',ll1); 
                    my_fprintf('Warp:\t%g\t%g\t%g\t%g\t%g :o)\t(%g)\n',ll1,lb,lh,llr1,llrb,armijo,verbose);
                    ll       = ll1;
                    llr      = llr1;
                    Twarp    = Twarp1;
                    def_done = 1;
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
            [x1,y1,z1] = defs(Twarp,zix,x0,y0,z0,M,buf(zix).code);
            tpm        = spm_sample_logpriors8(logtpm,x1,y1,z1);
            clear x1 y1 z1

            K1 = floor(sqrt(Kb));
            K2 = ceil(Kb/K1); 
            set(0,'CurrentFigure',fig{2});                                        
            for i=1:Kb     
                subplot(K1,K2,i);
                slice(buf(zix).code>0) = tpm{i};
                imagesc(slice'); axis image xy off; colormap(gray);
                title(['TPM, k=' num2str(i)]);
            end 
            clear tpm
            drawnow
        end
    end

    if ~(abs(ll-ooll)>2*tol1*nm) && iter>9
        % Finished
        break
    end
end

if dotpm       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For estimating a template using maximum likelihood
    %------------------------------------------------------------
    
    % Compute likelihoods
    %-----------------------------------------------------------
    px = zeros([prod(d(1:2)) d(3) Kb],'single');
    for z=1:length(z0)
        if ~buf(z).Nm, continue; end
              
        qt     = log_likelihoods(K,buf(z),[],mog);                
        max_qt = max(qt,[],2);
        msk    = buf(z).code>0;
        q      = zeros([buf(z).Nm Kb]);
        for k1=1:Kb
            for k=find(lkp==k1)
                q(:,k1) = q(:,k1) + exp(qt(msk,k) - max_qt(msk));
            end
            px(msk,z,k1) = single(q(:,k1));
        end
        clear qt max_qt q msk
    end          
    px = reshape(px,[d Kb]);    
    clear buf

    if ~isempty(fig{4}) 
        % Display likelihoods in subject space
        set(0,'CurrentFigure',fig{4});      
        for i=1:Kb
            subplot(2,Kb,i);
            imagesc(px(:,:,zix,i)'); axis image xy off; colormap(gray);
            title(['px, k=' num2str(i)]);
        end  
    end
    
    % Push likelihoods (in subject space) to template space 
    %-----------------------------------------------------------
    t      = identity(d) + bsxfun(@times,Twarp,1./sk4); 
    MM     = M*MT;
    t      = affine_transf(MM,t);     
    if dtpm(3)==1, t(:,:,:,3) = 1; end
    [px,c] = spm_diffeo('push',px,t,dtpm);
    clear t                         
    
    % Mask
    msk = isfinite(px(:,:,:,1));    
    msk = bsxfun(@and,msk,c>0);
    clear c
    
    % Calculate a bounding box from mask (TODO: implement own version of imBoundingBox)    
    bb = imBoundingBox(msk); % [y x z]    
    clear msk
    
    % Get bounding box info
    [dm,ix_x,ix_y,ix_z] = bb_info(bb,dtpm);
    ntpm                = prod(dm);
    
    % 'Crop' likelihoods
    px = px(ix_x,ix_y,ix_z,:);   
    
%     % Renormalise
%     px = max(px,eps);
%     px = bsxfun(@rdivide,px,sum(px,4));    
    
    if ~isempty(fig{4}) 
        % Display likelihoods pushed to template space
        set(0,'CurrentFigure',fig{4});      
        for i=(Kb + 1):2*Kb
            subplot(2,Kb,i);
            imagesc(px(:,:,floor(dm(3)/2 + 1),i - Kb)'); axis image xy off; colormap(gray);
            title(['push(px), k=' num2str(i - Kb)]);
        end  
    end           
    
    % Solve dL/dtpm=0 for tpm (where L is the objective function and tpm is the template)  
    %-----------------------------------------------------------
       
    % log(tpm) [Kb ntpm]
    logtpm_bb = zeros([dm Kb],'single'); 
    for k=1:Kb
        logtpm_bb(:,:,:,k) = logtpm.dat{k}(ix_x,ix_y,ix_z);
    end
    clear logtpm        
        
    % exp(log(wp) + log(tpm)) = wp*tpm
    logtpm_bb = reshape(logtpm_bb,[ntpm Kb])';
    logwp     = repmat(log(wp)',1,ntpm);    
    wtpm      = bsxfun(@plus,logtpm_bb,logwp);    
    wtpm      = exp(wtpm);
    clear logtpm_bb

    % wp*tpm*P(x|theta)    
    px     = reshape(px,[ntpm Kb])'; 
    wtpmpx = bsxfun(@times,wtpm,px);
    clear px

    % wp*tpm*P(x|theta)/sum(wp*tpm*P(x|theta),k)
    tpmnum = bsxfun(@rdivide,wtpmpx,sum(wtpmpx,1)); 
    tpmnum(~isfinite(tpmnum)) = 0;
    clear wtpmpx                 
    
    % Save numerator of update
    if exist(obj.pth_tpmnum,'file')==2
        delete(obj.pth_tpmnum);
    end
    create_nii(obj.pth_tpmnum,zeros([dm Kb],'single'),eye(4),'float32','tpmnum')                    
    
    Nii              = nifti(obj.pth_tpmnum);
    tpmnum           = reshape(tpmnum',[dm Kb]); 
    Nii.dat(:,:,:,:) = tpmnum;
    clear Nii tpmnum
    
    % wp/sum(wp*tpm,k)
    wtpm   = 1./sum(wtpm,1);     
    tpmden = bsxfun(@times,wtpm,exp(logwp));
    tpmden(~isfinite(tpmden)) = 0; 
    clear wtpm
    
    % Save denominator of update
    if exist(obj.pth_tpmden,'file')==2
        delete(obj.pth_tpmden);
    end
    create_nii(obj.pth_tpmden,zeros([dm Kb],'single'),eye(4),'float32','tpmden') 
    
    Nii              = nifti(obj.pth_tpmden);
    tpmden           = reshape(tpmden',[dm Kb]); 
    Nii.dat(:,:,:,:) = tpmden;
    clear Nii tpmden        
end

% Save the results
Nii              = nifti(obj.pth_def);
Nii.dat(:,:,:,:) = Twarp; 
clear Nii

obj.def_done = def_done;
obj.Tbias    = {chan(:).T};
obj.mg       = mg;
obj.wp       = wp;
obj.ll       = ll;
obj.nm       = nm;
obj.mog      = mog;
obj.bb       = bb;
    
return;
%==========================================================================
      
%==========================================================================
function [x1,y1,z1] = defs(Twarp,z,x0,y0,z0,M,code)
x1a = x0    + double(Twarp(:,:,z,1));
y1a = y0    + double(Twarp(:,:,z,2));
z1a = z0(z) + double(Twarp(:,:,z,3));
if nargin>=7  
    msk = code>0;
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end
[x1,y1,z1] = affine_transf(M,x1a,y1a,z1a);
if numel(z0)==1, z1 = ones(size(z1)); end
return;
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%==========================================================================
                 
%==========================================================================
function [Kb,K,lkp,mg,mog,mom] = split_mog(lkp,mog,omom)
% A crude heuristic to replace a single Gaussian by a bunch of Gaussians
% If there is only one Gaussian, then it should be the same as the
% original distribution.
vb  = mog.vb;
K   = numel(lkp);
Kb  = max(lkp);
mg  = ones(1,K)/K;

if vb           
    m = mog.po.m;
    b = mog.po.b;
    W = mog.po.W;
    n = mog.po.n;   

    m0 = mog.pr.m;
    b0 = mog.pr.b;
    W0 = mog.pr.W;
    n0 = mog.pr.n;
    
    N = size(m,1);
    
    po.m = zeros(N,K);
    po.b = zeros(1,K);
    po.W = zeros(N,N,K);
    po.n = zeros(1,K);

    pr.m = zeros(N,K);
    pr.b = zeros(1,K);
    pr.W = zeros(N,N,K);
    pr.n = zeros(1,K);
else    
    mn0 = mog.mn;
    vr0 = mog.vr;        
    
    N = size(mn0,1);
    
    mn  = zeros(N,K);                    
    vr  = zeros(N,N,K);        
end

omom = mom_John2Bishop(omom);
mom  = mom_struct(K,N); 
for i=2:2^N           
    for k=1:Kb
        if omom(i).s0(k)>eps*eps
            kk = sum(lkp==k);
            w  = 1./(1 + exp(-(kk-1)*0.25)) - 0.5; 

            s0 = repmat(omom(i).s0(k)/kk,1,kk);        
            s1 = sqrtm(omom(i).S2(:,:,k))*randn(N,kk)*w + repmat(omom(i).s1(:,k),[1,kk]);
            S2 = repmat(omom(i).S2(:,:,k)*(1 - w),[1,1,kk]);                                      
            
            mom(i).s0(lkp==k)     = s0;
            mom(i).s1(:,lkp==k)   = s1;
            mom(i).S2(:,:,lkp==k) = S2;
        end
    end  
end
mom = mom_Bishop2John(mom);

for k=1:Kb
    kk 	       = sum(lkp==k);
    mg(lkp==k) = 1/kk;
    w  		   = 1./(1 + exp(-(kk-1)*0.25)) - 0.5; 

    if vb
        pr.b(lkp==k)     = b0(k);
        pr.m(:,lkp==k)   = repmat(m0(:,k),[1,kk]);
        pr.n(lkp==k)     = n0(k);
        pr.W(:,:,lkp==k) = repmat(W0(:,:,k),[1 1 kk]); 

        vr0              = inv(n0(k)*W0(:,:,k));
        pr0              = inv(vr0*(1 - w));                        
        pr.b(lkp==k)     = b0(k)/kk;
        pr.m(:,lkp==k)   = sqrtm(vr0)*randn(N,kk)*w + repmat(m0(:,k),[1,kk]);
        pr.n(lkp==k)     = n0(k)/kk;
        pr.W(:,:,lkp==k) = repmat(pr0/n0(k),[1 1 kk]); 
    else
        mn(:,lkp==k)   = sqrtm(vr0(:,:,k))*randn(N,kk)*w + repmat(mn0(:,k),[1,kk]);
        vr(:,:,lkp==k) = repmat(vr0(:,:,k)*(1 - w),[1,1,kk]);
    end                                        
end  

if vb
    mog.pr = pr; 
    mog.po = po;
else
    mog.mn = mn; 
    mog.vr = vr;
end
%==========================================================================

%==========================================================================
function mom = mom_from_buf(buf,Kb,N,z0,uniform)
if uniform
    mom = spm_kmeans2mom(buf,Kb);
else
    mom = mom_struct(Kb,N); 
    for z=1:length(z0)     
        if ~buf(z).Nm, continue; end

        cr = zeros(numel(buf(z).code),N);
        for n=1:N
            cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n});
        end

        % Responsibilities from TPMs
        q = zeros(numel(buf(z).code),Kb);
        for k=1:Kb
            q(buf(z).code>0,k) = double(buf(z).dat(:,k));
        end

        % Calculate sufficient statistics                
        mom = spm_SuffStats(cr,q,buf(z).code,mom);                                                          
    end   
end 
%==========================================================================

%==========================================================================
function count = my_fprintf(varargin)
if varargin{end}
    count = fprintf(varargin{1:end - 1});
else
    count = 0;
end
%==========================================================================