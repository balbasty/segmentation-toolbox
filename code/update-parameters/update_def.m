function [ll,llr,buf,Twarp,L,armijo] = update_def(ll,llrb,llr,buf,mg,gmm,wp,lkp,Twarp,sk4,M,MT,tpm,x0,y0,z0,param,iter,fig,L,print_ll,tot_S,armijo,wp_lab)
nz   = numel(buf);
N    = numel(buf(1).f);
K    = numel(lkp.part);
Kb   = max(lkp.part);
d    = size(buf(1).msk{1});

Alpha = zeros([size(x0),nz,6],'single');
Beta  = zeros([size(x0),nz,3],'single');
for z=1:nz
    if ~buf(z).Nm, continue; end

    % Deformations from parameters
    msk1       = buf(z).code>0;
    [x1,y1,z1] = make_deformation(Twarp,z,x0,y0,z0,M,msk1);

    % Tissue probability map and spatial derivatives
    [b,db1,db2,db3] = spm_sample_logpriors(tpm,x1,y1,z1);
    clear x1 y1 z1

    % Adjust for tissue weights
    s   = zeros(size(b{1}));
    ds1 = zeros(size(b{1}));
    ds2 = zeros(size(b{1}));
    ds3 = zeros(size(b{1}));
    for k1=1:Kb
        b{k1}   = wp(k1)*b{k1};
        db1{k1} = wp(k1)*db1{k1};
        db2{k1} = wp(k1)*db2{k1};
        db3{k1} = wp(k1)*db3{k1};
        s       =  s  + b{k1};
        ds1     = ds1 + db1{k1};
        ds2     = ds2 + db2{k1};
        ds3     = ds3 + db3{k1};
    end
    for k1=1:Kb
        b{k1}   = b{k1}./s;
        db1{k1} = (db1{k1}-b{k1}.*ds1)./s;
        db2{k1} = (db2{k1}-b{k1}.*ds2)./s;
        db3{k1} = (db3{k1}-b{k1}.*ds3)./s;
    end
    clear s ds1 ds2 ds3

    % Rotate gradients (according to initial affine registration) and
    % compute the sums of the tpm and its gradients, times the likelihoods
    % (from buf.dat).
    p   = zeros(buf(z).Nm,1) + eps;
    dp1 = zeros(buf(z).Nm,1);
    dp2 = zeros(buf(z).Nm,1);
    dp3 = zeros(buf(z).Nm,1);
    MM  = M*MT; % Map from sampled voxels to atlas data
    
    % Compute responsibilties        
    qt = latent(buf(z).f,buf(z).bf,mg,gmm,double(buf(z).dat),lkp,wp,buf(z).msk,buf(z).code,buf(z).labels,wp_lab);        
    q  = zeros(buf(z).Nm,Kb,'single');
    for k1=1:Kb
        for k=find(lkp.part==k1)
            q(:,k1) = q(:,k1) + qt(msk1,k);
        end
    end
    clear qt
    
    for k1=1:Kb
        pp  = double(q(:,k1));        
        p   = p   + pp.*b{k1};
        dp1 = dp1 + pp.*(MM(1,1)*db1{k1} + MM(2,1)*db2{k1} + MM(3,1)*db3{k1});
        dp2 = dp2 + pp.*(MM(1,2)*db1{k1} + MM(2,2)*db2{k1} + MM(3,2)*db3{k1});
        dp3 = dp3 + pp.*(MM(1,3)*db1{k1} + MM(2,3)*db2{k1} + MM(3,3)*db3{k1});
    end
    clear b db1 db2 db3

    % Compute first and second derivatives of the matching term.  Note that
    % these can be represented by a vector and tensor field respectively.
    tmp       = zeros(d(1:2));
    tmp(msk1) = dp1./p; dp1 = tmp; dp1 = reshape(dp1,d(1:2));
    tmp(msk1) = dp2./p; dp2 = tmp; dp2 = reshape(dp2,d(1:2));
    tmp(msk1) = dp3./p; dp3 = tmp; dp3 = reshape(dp3,d(1:2));

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
    
    % Adjust for label weight
    Beta  = (1 - wp_lab)*Beta;
    Alpha = (1 - wp_lab)*Alpha;
end

if tot_S==1
    % Heavy-to-light regularisation
    scal     = 2^max(10 - iter,0);       
    param(6) = param(6)*scal^2;
end

% Add in the first derivatives of the prior term
Beta = Beta  + spm_diffeo('vel2mom',bsxfun(@times,Twarp,1./sk4),param);

% Gauss-Newton increment
Update = bsxfun(@times,spm_diffeo('fmg',Alpha,Beta,[param 2 2]),sk4);
clear Alpha Beta

% Line search to ensure objective function improves
for line_search=1:12
    Twarp1 = Twarp - armijo*Update; % Backtrack if necessary
    if nz==1
        Twarp1(:,:,:,end) = single(0);
    end
    
    % Recompute objective function
    llr1 = -0.5*sum(sum(sum(sum(Twarp1.*bsxfun(@times,spm_diffeo('vel2mom',bsxfun(@times,Twarp1,1./sk4),param),1./sk4)))));
    ll1  = llr1 + llrb;

    mom = moments_struct(K,N);
    for z=1:nz
        if ~buf(z).Nm, continue; end

        msk1                          = buf(z).code>0;
        [x1,y1,z1]                    = make_deformation(Twarp1,z,x0,y0,z0,M,msk1);
        b                             = spm_sample_logpriors(tpm,x1,y1,z1);
        for k1=1:Kb, buf(z).dat(:,k1) = b{k1}; end
        clear x1 y1 z1

        cr                             = NaN(numel(buf(z).msk),N);
        for n=1:N, cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end 

        [q,dll] = latent(buf(z).f,buf(z).bf,mg,gmm,double(buf(z).dat),lkp,wp,buf(z).msk,buf(z).code,buf(z).labels,wp_lab,cr);
        ll1     = ll1 + dll;

        mom = spm_SuffStats(cr,q,mom,buf(z).code);
        clear q cr b
    end           

    % Compute missing data and VB components of ll
    dll = spm_VBGaussiansFromSuffStats(mom,gmm);
    ll1 = ll1 + sum(sum(dll)); 

    if ll1<ll
        % Still not better, so keep searching inwards.
        my_fprintf('Warp:\t%g\t%g\t%g :o(\t(%g)\n', ll1, llr1,llrb,armijo,print_ll);
        
        armijo = armijo*0.75;
        
        if line_search==12
            L{1}(end + 1) = ll;
            L{3}(end + 1) = llr;
            
            for z=1:nz
                % Revert to previous deformation
                if ~buf(z).Nm, continue; end
                
                msk1       = buf(z).code>0;
                [x1,y1,z1] = make_deformation(Twarp,z,x0,y0,z0,M,msk1);
                b          = spm_sample_logpriors(tpm,x1,y1,z1);
                clear x1 y1 z1
                
                for k1=1:Kb, buf(z).dat(:,k1) = b{k1}; end                
            end
        end   
    else
        % Better.  Accept the new solution.            
        my_fprintf('Warp:\t%g\t%g\t%g :o)\t(%g)\n', ll1, llr1,llrb,armijo,print_ll);
        
        armijo = min(armijo*1.2,1);
        
        ll    = ll1;
        llr   = llr1;
        Twarp = Twarp1;
        
        L{1}(end + 1) = ll;
        L{3}(end + 1) = llr;
                
        debug_view('convergence',fig{4},lkp,buf,L);
        
        break
    end
end
%========================================================================== 

%==========================================================================
function [x1,y1,z1] = make_deformation(Twarp,z,x0,y0,z0,M,msk)
x1a = x0    + double(Twarp(:,:,z,1));
y1a = y0    + double(Twarp(:,:,z,2));
z1a = z0(z) + double(Twarp(:,:,z,3));
if nargin>=7
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
if numel(z0)==1
   z1 = ones(size(z1));
end
return;
%==========================================================================