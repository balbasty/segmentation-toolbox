function mom = spm_kmeans2mom(buf,K,verbose)
% Generate initial estimates for the parameters of a Gaussian mixture model
% using the k-means algorithm
%
% FORMAT mom = spm_kmeans2mom(buf,K,verbose)
%     
if nargin<3, verbose = 0; end

N = numel(buf(1).f);
d = [size(buf(1).msk{1}) numel(buf)];

F = NaN([prod(d(1:2)) d(3) N],'single');
for z=1:numel(buf)
    for n=1:N
        F(buf(z).msk{n},z,n) = buf(z).f{n};
    end
end
F = reshape(F,[d N]);    

if verbose
    % Display input image(s)    
   figure(666);
   for n=1:N
      subplot(1,N,n);
      imagesc(F(:,:,floor(d(3)/2) + 1,n)'); axis image xy off; colormap(gray);
   end
   drawnow
end

% Label images using k-means
F = reshape(F,[prod(d) N]);
Q = label_data(F,K,d);

if verbose
    % Display estimated labels    
    figure(667);
    tmp = reshape(Q,[d K]);                   
    K1  = floor(sqrt(K));
    K2  = ceil(K/K1); 
    for k=1:K
      subplot(K1,K2,k);
      imagesc(tmp(:,:,floor(d(3)/2) + 1,k)'); axis image xy off; colormap(gray);
    end
    clear tmp
    drawnow
end

% Generate estimates of MoG parameters
F   = reshape(F,[d N]);
Q   = reshape(Q,[d K]);
mom = mom_struct(K,N);  
for z=1:d(3)    
    f   = reshape(double(F(:,:,z,:)),[prod(d(1:2)) N]);
    q   = reshape(double(Q(:,:,z,:)),[prod(d(1:2)) K]);
    mom = spm_SuffStats(f,q,buf(z).code,mom);
end
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,d)
% w = warning('query','last')
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',1000);

labels = kmeans(f,single(K),...
                'Distance','sqeuclidean',...
                'Start','plus',...
                'Replicates',5,...
                'Options',opts);
labels = single(labels);

labels                    = labels';
labels(~isfinite(labels)) = K + 1;

nlabels = zeros([prod(d), K + 1],'single');

idx          = sub2ind(size(nlabels),1:prod(d),labels);
nlabels(idx) = 1;
clear labels

idx = nlabels(:,K + 1) == 1;    
for k=1:K
   nlabels(idx,k) = 1/K;
end
nlabels(:,K + 1) = [];

warning('on','stats:kmeans:MissingDataRemoved')
warning('on','stats:kmeans:FailedToConvergeRep')
%==========================================================================