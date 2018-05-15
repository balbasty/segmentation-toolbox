function mom = kmeans2mom(buf,lkp,mn,mx,obj)
N  = numel(buf(1).f); 
Kb = max(lkp.part); 

C = spm_misc('linspace_vec',mn,mx,Kb);
if size(C,1)~=Kb, C = C'; end    

F = [];
for z=1:numel(buf)          
    if ~numel(buf(z).f{1}), continue; end 
    
    cr                 = zeros([numel(buf(z).f{1}) N],'single');
    for n=1:N, cr(:,n) = buf(z).f{n}.*buf(z).bf{n}; end
    F                  = [F;cr];
    clear cr    
end

Q = label_data(F,Kb,C,obj.segment.kmeans_dist);

mom = moments_struct(Kb,N);  
for k=1:Kb
    mom(end).s0(1,k)   = mom(end).s0(1,k)   + sum(Q(:,k));
    mom(end).s1(:,k)   = mom(end).s1(:,k)   + F(:,:)'*Q(:,k);
    mom(end).S2(:,:,k) = mom(end).S2(:,:,k) + bsxfun(@times,Q(:,k),F(:,:))'*F(:,:);
end
% mom(end).s1./mom(end).s0

if ~isempty(obj.segment.init_clust)
    if strcmp(obj.segment.init_clust,'random')
        % Randomise sufficient statistics        
        ix = randperm(Kb);        
    elseif  strcmp(obj.segment.init_clust,'total')
        % Sort sufficient statistics according total amount of each class (0th moment)
        [~,ix] = sort(mom(end).s0);  
    elseif  strcmp(obj.segment.init_clust,'magnitude')
        % Reverse sufficient statistics
        [~,ix] = sort(sum(sqrt(mom(end).s1.^2),1));
    elseif  strcmp(obj.segment.init_clust,'reverse')
        % Reverse sufficient statistics
        ix = fliplr(1:Kb);   
    elseif  strcmp(obj.segment.init_clust,'user')
        % Used-defined
        ix = obj.segment.class_ix;        
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