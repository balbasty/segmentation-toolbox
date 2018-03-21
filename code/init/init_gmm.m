function [gmm,buf] = init_gmm(obj,N,buf,vr0,mn,mx)
gmm            = obj.gmm;
gmm.init_clust = obj.init_clust;
gmm.ml         = obj.do_ml;

gmm.vr0 = vr0;
gmm.min = mn;
gmm.max = mx;

if (~isfield(gmm,'mn') && ~isfield(gmm,'vr')) && ~isfield(gmm,'po')
    
    Kb  = max(obj.lkp.part);
    lkp = obj.lkp;
    
    % Compute moments
    %----------------------------------------------------------------------
    if obj.uniform
        % Uniform template provided, use the k-means algorithm to comppute
        % moments
        if isempty(lkp.lab)            
            mom = kmeans2mom(buf,numel(lkp.keep),mn,mx,obj.init_clust,obj.kmeans_dist,obj.kmeans_ix);        
        else                     
            % Extract non-labelled voxels
            msk  = ismember(lkp.keep,lkp.lab);
            lkp1 = lkp.keep(~msk);
            
            N    = numel(buf(1).f);
            buf1 = struct;
            for z=1:numel(buf)
                msk = sum(buf(z).labels,2)==0;
                buf1(z).img = cell(1,N);
                for n=1:N
                    buf1(z).img{n} = buf(z).img{n}(msk);
                end
            end
            
            mom1 = kmeans2mom(buf1,nnz(lkp1),mn,mx,obj.init_clust,obj.kmeans_dist);   
            clear buf1
            
            % Extract labelled voxels                                
            lkp2 = lkp.lab;
            
            mom2 = moments_struct(nnz(lkp2),N);
            for z=1:numel(buf)
                msk = sum(buf(z).labels,2)>0;
                
                cr          = zeros(nnz(msk),N);
                for n=1:N, 
                    cr(:,n) = double(buf(z).img{n}(msk)); 
                end  
                
                q          = zeros(nnz(msk),nnz(lkp2));
                for k=1:nnz(lkp2)
                    q(:,k) = double(buf(z).labels(msk,lkp.lab==k));
                end
                
                for k=1:nnz(lkp2)
                    mom2(end).s0(1,k)   = mom2(end).s0(1,k)   + sum(q(:,k));
                    mom2(end).s1(:,k)   = mom2(end).s1(:,k)   + cr(:,:)'*q(:,k);
                    mom2(end).S2(:,:,k) = mom2(end).S2(:,:,k) + bsxfun(@times,q(:,k),cr(:,:))'*cr(:,:);
                end
            end
            
            % Combine moments
            mom = moments_struct(numel(lkp.keep),N);
            cnt = 1;
            for k=1:numel(lkp.keep)
                if lkp2(k)
                    mom(end).s0(1,k)   = mom2(end).s0(1,lkp2(k));
                    mom(end).s1(:,k)   = mom2(end).s1(:,lkp2(k));
                    mom(end).S2(:,:,k) = mom2(end).S2(:,:,lkp2(k));                    
                else
                    mom(end).s0(1,k)   = mom1(end).s0(1,cnt);
                    mom(end).s1(:,k)   = mom1(end).s1(:,cnt);
                    mom(end).S2(:,:,k) = mom1(end).S2(:,:,cnt);
                    cnt                = cnt + 1;
                end
            end
        end
        buf = rmfield(buf,'img');
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
