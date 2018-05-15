function obj = spm_segment(obj,fig)

% Parameters, etc.
%-----------------------------------------------------------------------
V         = obj.image; 
N         = numel(V);
d0        = V(1).dim(1:3);
vx        = spm_misc('vxsize',V(1).mat);
sk        = max([1 1 1],round(obj.segment.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
d         = [size(x0) length(z0)];
lkp       = obj.segment.lkp;
Kb        = max(lkp.part);
tol1      = obj.segment.tol1;
wp_reg    = obj.segment.wp_reg;
wp_l      = obj.segment.wp_l;
modality  = obj.modality;
niter     = obj.segment.niter;
nitgmm    = obj.segment.nitgmm;
nsubit    = obj.segment.nsubit;
do_bf     = obj.segment.do_bf;
if obj.uniform, obj.segment.do_def = false; end
do_def    = obj.segment.do_def;
do_wp     = obj.segment.do_wp;
print_ll  = obj.segment.print_ll;
pth_vel   = obj.pth_vel;
Affine    = obj.Affine;
iter_template = obj.iter;
do_template = obj.do_template;
do_mg       = obj.segment.do_mg;

% Initialise weights
%-----------------------------------------------------------------------
if isfield(obj.segment,'wp') 
    wp = obj.segment.wp;
else
    wp = ones(1,Kb)/Kb;
end
wp(lkp.rem) = eps;
wp          = wp/sum(wp);

if isfield(obj.segment,'mg')
    mg = obj.segment.mg;   
else              
    mg       = ones(1,Kb);
end

% Load template
%-----------------------------------------------------------------------
tpm = spm_load_logpriors(obj.pth_template,wp);
M   = tpm.M\Affine*V(1).mat; % Affine matrix

% Fudge Factor - to (approximately) account for non-independence of voxels.
%-----------------------------------------------------------------------
ff = compute_fudge_factor(obj.fwhm,vx,sk);

% Load data into buffer
%-----------------------------------------------------------------------
[buf,nm,vr0,mn,mx,scl_bf] = init_buf(N,obj,V,x0,y0,z0,o,M,tpm);

% Initialise bias field
%-----------------------------------------------------------------------
[buf,chan,llrb] = init_bf(buf,obj,V,x0,y0,z0,ff,scl_bf);

% Initialise deformation and template
%-----------------------------------------------------------------------
[buf,param,MT,sk4,Twarp,llr] = init_def_and_dat(buf,obj,sk,vx,ff,d,fig,wp,x0,y0,z0,tpm,M);

% Initialise GMM
%-----------------------------------------------------------------------
[gmm,buf] = init_gmm(obj,buf,vr0,mn,mx);                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iterating
%------------------------------------------------------------
armijo = ones(1,2);
ll     = -Inf;
L      = {ll,llrb,llr};
for iter=1:niter

    if do_bf
        for subit=1:nsubit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate bias
            % Note that for multi-spectral data, the covariances among
            % channels are not computed as part of the second derivatives.
            % The aim is to save memory, and maybe make the computations
            % faster.
            %------------------------------------------------------------
                        
            % Estimate cluster parameters            
            [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_l,iter_template,do_mg);

            debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,wp_l);

            if subit>1 && ~((ll-oll)>2*tol1*nm), 
                break; 
            end
            oll = ll;
            
            [ll,llrb,buf,chan,L,armijo(1)] = update_bf(ll,llrb,llr,buf,mg,gmm,wp,lkp,chan,fig,L,print_ll,armijo(1),wp_l);                                          
        end
    end
    
    debug_view('bf',fig{2},lkp,buf,modality);
    
    if do_def
        for subit=1:nsubit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Estimate deformations
            %------------------------------------------------------------            
                        
            % Estimate cluster parameters            
            [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_l,iter_template,do_mg);

            debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,wp_l);

            if subit>1 && ~((ll-oll)>2*tol1*nm), 
                break; 
            end
            oll = ll;
        
            
            [ll,llr,buf,Twarp,L,armijo(2)] = update_def(ll,llrb,llr,buf,mg,gmm,wp,lkp,Twarp,sk4,M,MT,tpm,x0,y0,z0,param,iter,fig,L,print_ll,do_template,armijo(2),wp_l);
        end    
    end

    debug_view('template',fig{3},lkp,buf,wp);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate cluster parameters
    %------------------------------------------------------------
    [ll,mg,gmm,wp,L] = update_gmm(ll,llr,llrb,buf,mg,gmm,wp,lkp,wp_reg,iter,tol1,nm,nitgmm,do_wp,fig,L,print_ll,wp_l,iter_template,do_mg);    
    
    debug_view('responsibilities',fig{1},lkp,buf,gmm,mg,wp,wp_l);
        
    if iter>=10 && ~((ll-ooll)>2*tol1*nm)
        % Finished
        break
    end
    ooll = ll;
end
clear tpm

% Compute final moments, used for updating posteriors outside of main loop
%--------------------------------------------------------------------------
mom = compute_moments(buf,lkp,mg,gmm,wp,wp_l);   
clear buf

% For setting the DC component of all the bias fields so that they
% average to 0 (used for global signal normalisation).
%--------------------------------------------------------------------------
dc = zeros(1,N);
b1 = zeros(1,N);
b2 = zeros(1,N);
b3 = zeros(1,N);
for n=1:N    
    dc(n) = chan(n).T(1,1,1); 
    b1(n) = chan(n).B1(1,1); 
    b2(n) = chan(n).B2(1,1); 
    b3(n) = chan(n).B3(1,1); 
end 

% Save the results
%--------------------------------------------------------------------------
Nii = nifti(pth_vel);
for i=1:3
    Nii.dat(:,:,:,i) = Twarp(:,:,:,i); 
end
clear Twarp Nii

obj.Affine         = Affine;

obj.segment.lkp    = lkp;
obj.segment.MT     = MT;
obj.segment.Tbias  = {chan(:).T};
obj.segment.wp     = wp;
obj.segment.mg     = mg;
obj.segment.mom    = mom;
obj.segment.gmm    = gmm;
obj.segment.ll     = ll;
obj.segment.ll_all = [obj.segment.ll_all ll];

obj.segment.bf.dc = dc;
obj.segment.bf.b1 = b1;
obj.segment.bf.b2 = b2;
obj.segment.bf.b3 = b3;

return;
%==========================================================================