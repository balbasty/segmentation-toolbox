function mom = spm_kmeans2mom(buf,K,modality,ix_subj,verbose)
% Generate initial estimates for the parameters of a Gaussian mixture model
% using the k-means algorithm
%
% FORMAT mom = spm_kmeans2mom(buf,K,verbose)
%     
if nargin<5, verbose = true; end

N = numel(buf(1).f);
d = [size(buf(1).msk) numel(buf)];


mom = moments_struct(K,N);  
for z=1:d(3)      
    if ~buf(z).nm, continue; end
    
    cr                          = zeros(numel(buf(z).msk),N);
    for n=1:N, cr(buf(z).msk,n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end

    rng(1);
    q = label_data(cr,K,d);
         
%     figure(666);
%     K1 = floor(sqrt(K));
%     K2 = ceil(K/K1); 
%     for k=1:K
%       subplot(K1,K2,k);
%       tmp = reshape(q(:,k),d(1:2));
%       imagesc(tmp'); axis image xy off; colormap(gray);
%     end
%     clear tmp
%     drawnow
    
    mom = spm_SuffStats(cr,q,mom,buf(z).code);
end

% % Read data
% %--------------------------------------------------------------------------
% F = NaN([prod(d(1:2)) d(3) N],'single');
% for z=1:numel(buf)    
%     if ~buf(z).nm, continue; end
%     
%     for n=1:N
%         if isfield(buf(z),'bf'), F(buf(z).msk,z,n) = buf(z).f{n}.* buf(z).bf{n};
%         else,                    F(buf(z).msk,z,n) = buf(z).f{n};
%         end
%     end
% end
% F = reshape(F,[d N]);    
% 
% if verbose
%     % Display input image(s)    
%    figure(666);
%    for n=1:N
%       subplot(1,N,n);
%       imagesc(F(:,:,floor(d(3)/2) + 1,n)'); axis image xy off; colormap(gray);
%    end
%    drawnow
% end

% % Label data using k-means
% %--------------------------------------------------------------------------
% if strcmp(modality,'CT'), rng(ix_subj); end
%     
% F = reshape(F,[prod(d) N]);
% Q = label_data(F,K,d);
% Q = reshape(Q,[d K]);
% clear F

% % Calculate sufficient statistics
% %--------------------------------------------------------------------------
% mom = moments_struct(K,N);  
% for z=1:d(3)      
%     if ~buf(z).nm, continue; end
%     
%     cr                          = zeros(numel(buf(z).msk),N);
%     for n=1:N, cr(buf(z).msk,n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end
% 
%     q   = reshape(double(Q(:,:,z,:)),[prod(d(1:2)) K]);   
%     mom = spm_SuffStats(cr,q,mom,buf(z).code);
% end

% if strcmp(modality,'MRI')
%     % Sort sufficient statistics according to average Euclidean distance
%     %--------------------------------------------------------------------------
%     s1 = 0;
%     for i=2:numel(mom) 
%         s1 = s1 + sqrt(sum(mom(i).s1.^2,1));
%     end
%     s1 = s1/(numel(mom) - 1);
% 
%     [~,ix] = sort(s1);  
% elseif strcmp(modality,'CT')
    % Randomise sufficient statistics
    %--------------------------------------------------------------------------
    rng(ix_subj);
    ix = randperm(K);
    rng(1);
% end

for i=2:numel(mom) 
    mom(i).s0 = mom(i).s0(ix);
    mom(i).s1 = mom(i).s1(:,ix);
    mom(i).S2 = mom(i).S2(:,:,ix);
end

% if strcmp(modality,'CT'), rng(1); end

% if verbose
%     % Display estimated labels    
%     figure(667);    
%     Q  = Q(:,:,:,ix);        
%     K1 = floor(sqrt(K));
%     K2 = ceil(K/K1); 
%     for k=1:K
%       subplot(K1,K2,k);
%       imagesc(Q(:,:,floor(d(3)/2) + 1,k)'); axis image xy off; colormap(gray);
%     end
%     clear tmp
%     drawnow
% end
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,d)
% w = warning('query','last')
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',1000);

labels = kmeans(f,single(K),...
                'Distance','sqeuclidean',... % sqeuclidean, cityblock
                'Start','sample',...
                'Replicates',5,...
                'Options',opts);
labels = single(labels);

nlabels = zeros([numel(labels) K]);
for k=1:K, nlabels(:,k) = labels==k; end

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