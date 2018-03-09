function mom = kmeans2mom(buf,K,mn,mx,init_clust,kmeans_dist)
N = numel(buf(1).img);

C = linspace_vec(mn,mx,K);
if size(C,1)~=K, C = C'; end

F = [];
for z=1:numel(buf)          
    if ~numel(buf(z).img{1}), continue; end 
    
    cr                 = zeros([numel(buf(z).img{1}) N],'single');
    for n=1:N, cr(:,n) = buf(z).img{n}; end
    F                  = [F;cr];
    clear cr    
end

Q = label_data(F,K,C,kmeans_dist);

mom = moments_struct(K,N);  
for k=1:K
    mom(end).s0(1,k)   = mom(end).s0(1,k)   + sum(Q(:,k));
    mom(end).s1(:,k)   = mom(end).s1(:,k)   + F(:,:)'*Q(:,k);
    mom(end).S2(:,:,k) = mom(end).S2(:,:,k) + bsxfun(@times,Q(:,k),F(:,:))'*F(:,:);
end

if ~isempty(init_clust)
    if strcmp(init_clust,'random')
        % Randomise sufficient statistics        
        ix = randperm(K);        
    elseif  strcmp(init_clust,'total')
        % Sort sufficient statistics according total amount of each class (0th moment)
        [~,ix] = sort(mom(end).s0);  
    elseif  strcmp(init_clust,'reverse')
        % Reverse sufficient statistics
        ix = fliplr(1:K);   
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
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',1000);

labels = kmeans(f,K,...
                'Distance',kmeans_dist,...
                'Start',C,...
                'Options',opts);

nlabels                 = zeros([numel(labels) K],'single');
for k=1:K, nlabels(:,k) = labels(:)==k; end

warning('on','stats:kmeans:MissingDataRemoved')
warning('on','stats:kmeans:FailedToConvergeRep')
%==========================================================================