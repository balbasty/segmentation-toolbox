function [gmm,buf] = init_gmm(obj,buf,vr0,mn,mx)
gmm            = obj.segment.gmm;
gmm.init_clust = obj.segment.init_clust;
gmm.vr0        = vr0;
gmm.min        = mn;
gmm.max        = mx;
lkp            = obj.segment.lkp;
part           = lkp.part;
N              = size(vr0,1);
modality       = obj.modality;
do_template    = obj.do_template;

if ~isfield(gmm,'po')    
    Kb = max(part);
    K  = numel(part);            
    
    if isfield(gmm,'pr') && strcmp(modality,'CT') && do_template
       gmm.po = gmm.pr; 
       return
    end
    
    % Compute moments
    %----------------------------------------------------------------------
    if obj.uniform
        % Uniform template provided, use the k-means algorithm to compute
        % moments       
        mom = kmeans2mom(buf,lkp,mn,mx,obj); 
    else       
        % Use template to compute moments
        mom = template2mom(buf,lkp);
    end
    
    if Kb~=K
       % Adjust moments for when using multiple Gaussians per tissue
       mom1 = spm_misc('mom_John2Bishop',mom);
       m    = mom1(end).s1;
       vr   = mom1(end).S2;       
       
       tmp_gmm.po.m = m;
       tmp_gmm.po.n = (N - 0.999)*ones(1,Kb);
       tmp_gmm.po.b = ones(1,Kb);
       tmp_gmm.po.W = zeros(N,N,Kb);
       for k=1:Kb
           tmp_gmm.po.W(:,:,k) = inv(vr(:,:,k))/tmp_gmm.po.n(k);
       end
       tmp_gmm.pr.m = tmp_gmm.po.m;
       tmp_gmm.pr.n = tmp_gmm.po.n;
       tmp_gmm.pr.b = tmp_gmm.po.b;
       tmp_gmm.pr.W = tmp_gmm.po.W;       
       
       tmp_gmm = more_gaussians(tmp_gmm,part);
       
       mom = moments_struct(K,N);
       for k=1:K
           mom(end).s0(k)     = mom1(end).s0(lkp.part(k))/sum(lkp.part==lkp.part(k));
           mom(end).s1(:,k)   = tmp_gmm.po.m(:,k);
           mom(end).S2(:,:,k) = inv(tmp_gmm.po.n(k)*tmp_gmm.po.W(:,:,k));
       end

       mom = spm_misc('mom_Bishop2John',mom);
    end
        
    if ~isfield(gmm,'pr')    
        % Create priors
        mom1 = spm_misc('mom_John2Bishop',mom);
    
        vr1 = zeros(N,N);
        for k=1:K
            vr1 = vr1 + (mom1(end).S2(:,:,k) - mom1(end).s1(:,k)*mom1(end).s1(:,k)'/mom1(end).s0(k)); 
        end
        vr1 = (vr1 + N*vr0)/(sum(mom1(end).s0) + N);  
        
        m0 = mom1(end).s1;
        b0 = ones(1,K);
        n0 = (N - 0.999)*ones(1,K);
        W0 = zeros(N,N,K);
        for k=1:K
            W0(:,:,k) = inv(vr1)/n0(k);
        end
        
        gmm.pr.m = m0;
        gmm.pr.b = b0;
        gmm.pr.n = n0;
        gmm.pr.W = W0;
        
        gmm.pr.part = part; % just so that the partition can be loaded from the prior object
    else
        if numel(gmm.pr.b)~=K
           error('numel(gmm.pr.b)~=K') 
        end
    end
    
    % Estimate GMM parameters
    %----------------------------------------------------------------------
    [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);
end
%==========================================================================    

%==========================================================================    
function mom = template2mom(buf,lkp)
Kb = max(lkp.part);
nz = numel(buf);
N  = numel(buf(1).f);
  
% Compute responsibilities from template
mom = moments_struct(Kb,N);
for z=1:nz
    if ~buf(z).Nm, continue; end
            
    % Get BX (bias-field x image)
    cr                             = NaN(numel(buf(z).msk{1}),N);
    for n=1:N, cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end    
        
    msk1 = buf(z).code>0;
    q    = NaN(numel(buf(z).msk{1}),Kb);        
    for k=1:Kb
        q(msk1,k) = double(buf(z).dat(:,k));
    end
        
    % Update sufficient statistics
    mom = spm_SuffStats(cr,q,mom,buf(z).code);
end
%==========================================================================    