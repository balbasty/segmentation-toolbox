function [gmm,buf] = init_gmm(obj,N,buf,vr0,mn,mx)
gmm            = obj.segment.gmm;
gmm.init_clust = obj.segment.init_clust;
gmm.ml         = obj.segment.do_ml;

gmm.vr0 = vr0;
gmm.min = mn;
gmm.max = mx;

if (~isfield(gmm,'mn') && ~isfield(gmm,'vr')) && ~isfield(gmm,'po')
    
    lkp = obj.segment.lkp;
    Kb  = max(lkp.part);
        
    % Compute moments
    %----------------------------------------------------------------------
    if obj.uniform
        % Uniform template provided, use the k-means algorithm to compute
        % moments       
        mom = kmeans2mom(buf,lkp,mn,mx,obj);        
        
        buf = rmfield(buf,'img');        
        buf = rmfield(buf,'labels_full');
    else       
        % Use template to compute moments
        mom = compute_moments(buf,Kb);
    end
    
    % Estimate GMM parameters
    %----------------------------------------------------------------------
    if gmm.ml
        vr1     = zeros(N,N);
        [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm,vr1);  
    else       
        [~,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);  
    end    
    
    % Adjust GMM parameters when forcing responsibilities for some classes to zero
    %----------------------------------------------------------------------    
    ogmm = gmm;
    gmm  = ogmm;
    kk   = 1;
    for k=1:Kb
        if any(lkp.rem==k)
            if gmm.ml
                gmm.mn(:,k)   = zeros(size(ogmm.mn(:,1)));
                gmm.vr(:,:,k) = eye(size(ogmm.vr(:,:,1)));
            else
                gmm.pr.n(k)     = (N - .999)*ones(size(ogmm.pr.n(1)));
                gmm.pr.b(k)     = ones(size(ogmm.pr.b(1)));
                gmm.pr.m(:,k)   = zeros(size(ogmm.pr.m(:,1)));
                gmm.pr.W(:,:,k) = eye(size(ogmm.pr.W(:,:,1)));
                
                gmm.po.n(k)     = (N - .999)*ones(size(ogmm.po.n(1)));
                gmm.po.b(k)     = ones(size(ogmm.po.b(1)));
                gmm.po.m(:,k)   = zeros(size(ogmm.po.m(:,1)));
                gmm.po.W(:,:,k) = eye(size(ogmm.po.W(:,:,1)));
            end
        else
            if gmm.ml
                gmm.mn(:,k)   = ogmm.mn(:,kk);
                gmm.vr(:,:,k) = ogmm.vr(:,:,kk);
            else
                gmm.pr.n(k)     = ogmm.pr.n(kk);
                gmm.pr.b(k)     = ogmm.pr.b(kk);
                gmm.pr.m(:,k)   = ogmm.pr.m(:,kk);
                gmm.pr.W(:,:,k) = ogmm.pr.W(:,:,kk);   
                
                gmm.po.n(k)     = ogmm.po.n(kk);
                gmm.po.b(k)     = ogmm.po.b(kk);
                gmm.po.m(:,k)   = ogmm.po.m(:,kk);
                gmm.po.W(:,:,k) = ogmm.po.W(:,:,kk);                
            end  
            kk = kk + 1;
        end
    end
end
%======================================================================= 
