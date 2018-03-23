function mom = kmeans2mom(buf,lkp,mn,mx,obj)
N = numel(buf(1).img);
K = numel(lkp.keep);

C = linspace_vec(mn,mx,K);
if size(C,1)~=K, C = C'; end
    
if ~isempty(lkp.lab)    
    Nii    = nifti(obj.labels.fname);
    labels = uint8(Nii.dat(:,:,:));        
    if min(labels(:))==0
        labels = labels + 1;
    end
    
    avg_int = zeros(K,N);
    
    for n=1:N
        Nii = nifti(obj.image(n).fname);
        img = single(Nii.dat(:,:,:));
    
        for k=lkp.lab
            if k==0, continue; end

            msk          = labels==k;
            avg_int(k,n) = mean(img(msk(:)));
        end
    end
    clear labels img
    
    for k=lkp.lab
        if k==0, continue; end
        
        C(k,:) = avg_int(k,:);
    end
    
    C(lkp.lab==0,:) = linspace_vec(mn,mx,nnz(lkp.lab==0))';
end

F = [];
for z=1:numel(buf)          
    if ~numel(buf(z).img{1}), continue; end 
    
    cr                 = zeros([numel(buf(z).img{1}) N],'single');
    for n=1:N, cr(:,n) = buf(z).img{n}; end
    F                  = [F;cr];
    clear cr    
end

Q = label_data(F,K,C,obj.segment.kmeans_dist);

mom = moments_struct(K,N);  
for k=1:K
    mom(end).s0(1,k)   = mom(end).s0(1,k)   + sum(Q(:,k));
    mom(end).s1(:,k)   = mom(end).s1(:,k)   + F(:,:)'*Q(:,k);
    mom(end).S2(:,:,k) = mom(end).S2(:,:,k) + bsxfun(@times,Q(:,k),F(:,:))'*F(:,:);
end
% mom(end).s1./mom(end).s0

if ~isempty(obj.segment.init_clust)
    if strcmp(obj.segment.init_clust,'random')
        % Randomise sufficient statistics        
        ix = randperm(K);        
    elseif  strcmp(obj.segment.init_clust,'total')
        % Sort sufficient statistics according total amount of each class (0th moment)
        [~,ix] = sort(mom(end).s0);  
    elseif  strcmp(obj.segment.init_clust,'reverse')
        % Reverse sufficient statistics
        ix = fliplr(1:K);   
    elseif  strcmp(obj.segment.init_clust,'user')
        % Used-defined
        ix = obj.segment.kmeans_ix;        
    end
    
    for i=1:numel(mom) 
        mom(i).s0 = mom(i).s0(ix);
        mom(i).s1 = mom(i).s1(:,ix);
        mom(i).S2 = mom(i).S2(:,:,ix);
    end
end
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,C,kmeans_dist)
% w = warning('query','last')
warning('off','stats:kmeans:FailedToConverge')

opts = statset('MaxIter',1000);

labels = kmeans(f,K,...
                'Distance',kmeans_dist,...
                'Start',C,...
                'Options',opts);

nlabels                 = zeros([numel(labels) K],'single');
for k=1:K, nlabels(:,k) = labels(:)==k; end

warning('on','stats:kmeans:FailedToConverge')
%==========================================================================