function [wp,mn,vr,Q] = spm_InitGaussians(buf,K,verbose)
% Generate initial estimates for the parameters of a Gaussian mixture model
% using the k-means algorithm
%
% FORMAT [wp,mn,vr,Q] = spm_InitGaussians(buf,K,verbose)
%     
if nargin<3, verbose = false; end

N = numel(buf(1).f);
d = [size(buf(1).msk) numel(buf)];

f = NaN([prod(d(1:2)) d(3) N],'single');
for z=1:numel(buf)
    for n=1:N
        f(buf(z).msk,z,n) = buf(z).f{n};
    end
end
f = reshape(f,[d N]);    

if verbose
    % Display input image(s)    
   figure(get_nbr_figs + 1);
   for n=1:N
      subplot(1,N,n);
      imagesc(f(:,:,floor(d(3)/2) + 1,n)'); axis image xy off; colormap(gray);
   end
   drawnow
end

% Label images using k-means
f = reshape(f,[prod(d) N]);
Q = label_data(f,K,d);

if verbose
    % Display estimated labels    
    figure(get_nbr_figs + 1);
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
mom        = spm_SuffStats(f,Q);
[wp,mn,vr] = spm_GaussiansFromSuffStats(mom);
wp         = wp/sum(wp); 

% Sort estimates according to mean values of first image
[~,ix] = sort(mn(1,:),2);

mn = mn(:,ix);
vr = vr(:,:,ix);
wp = wp(:,ix);   

% Q  = uint8(255*Q(:,ix));
% Q  = reshape(Q,[d K]);
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,d)
% w = warning('query','last')
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',500);

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
   nlabels(idx,k) = NaN;
end
nlabels(:,K + 1) = [];

warning('on','stats:kmeans:MissingDataRemoved')
warning('on','stats:kmeans:FailedToConvergeRep')
%==========================================================================