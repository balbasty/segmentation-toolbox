function mom = kmeans2mom(buf,K,mn,mx,init_clust,kmeans_dist)
N = numel(buf(1).f);
d = [size(buf(1).msk{1}) numel(buf)];

C = linspace_vec(mn,mx,K);
if size(C,1)~=K, C = C'; end

mom = moments_struct(K,N);  
for z=1:d(3)          
    if sum(buf(z).code==(2^N - 1))<K, continue; end 
    
    cr                             = NaN(numel(buf(z).msk{1}),N);
    for n=1:N, cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end

    q = label_data(cr,K,d,C,kmeans_dist);

%     figure(666);
%     K1 = floor(sqrt(K));
%     K2 = ceil(K/K1); 
%     for k=1:K
%         subplot(K1,K2,k);
%         tmp = reshape(q(:,k),d(1:2));
%         imagesc(tmp'); axis image xy off; colormap(gray);
%     end
%     clear tmp
%     drawnow

    msk1 = buf(z).code==2^N - 1;
    for k=1:K
        mom(end).s0(1,k)   = mom(end).s0(1,k)   + sum(q(msk1,k));
        mom(end).s1(:,k)   = mom(end).s1(:,k)   + cr(msk1,:)'*q(msk1,k);
        mom(end).S2(:,:,k) = mom(end).S2(:,:,k) + bsxfun(@times,q(msk1,k),cr(msk1,:))'*cr(msk1,:);
    end
end

if strcmp(init_clust,'mean')
    % Sort sufficient statistics according mean (1st moment)
    %--------------------------------------------------------------------------
    s1 = sqrt(sum(mom(end).s1.^2,1));

    [~,ix] = sort(s1);  
    
elseif strcmp(init_clust,'total')    
    % Sort sufficient statistics according total amount of each class (0th moment)
    %--------------------------------------------------------------------------    
    [~,ix] = sort(mom(end).s0);  
    
elseif strcmp(init_clust,'random')
    % Randomise sufficient statistics
    %--------------------------------------------------------------------------    
    ix = randperm(K);    
else
    error('Wrong init!')
end

for i=1:numel(mom) 
    mom(i).s0 = mom(i).s0(ix);
    mom(i).s1 = mom(i).s1(:,ix);
    mom(i).S2 = mom(i).S2(:,:,ix);
end
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,d,C,kmeans_dist)
% w = warning('query','last')
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',1000);

labels = kmeans(f,K,...
                'Distance',kmeans_dist,...
                'Start',C,...
                'Options',opts);
msk    = isfinite(labels);

nlabels                   = NaN([numel(labels) K]);
for k=1:K, nlabels(msk,k) = labels(msk)==k; end

% labels                    = labels';
% labels(~isfinite(labels)) = K + 1;
% 
% nlabels = zeros([prod(d(1:2)), K + 1],'single');
% 
% idx          = sub2ind(size(nlabels),1:prod(d(1:2)),labels);
% nlabels(idx) = 1;
% clear labels
% 
% idx = nlabels(:,K + 1) == 1;    
% for k=1:K
%    nlabels(idx,k) = 0;
% end
% nlabels(:,K + 1) = [];

warning('on','stats:kmeans:MissingDataRemoved')
warning('on','stats:kmeans:FailedToConvergeRep')
%==========================================================================