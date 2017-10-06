function [nm,wp,mn,vr] = spm_InitGaussians(V,K,samp)
% Generate initial estimates for the parameters of a Gaussian mixture model
% using the k-means algorithm
%
% FORMAT [nm,wp,mn,vr] = spm_InitGaussians(V,K,samp)
%     
if nargin<3, samp=2; end

N = numel(V);

% Sub-sample images
d0        = V(1).dim(1:3);
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
d         = [size(x0) length(z0)];

f   = zeros([d N],'single');
msk = zeros(size(f),'logical');
for z=1:length(z0)
    for n=1:N
        f(:,:,z,n)   = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        msk(:,:,z,n) = get_msk(f(:,:,z,n));
    end
end
f   = reshape(f,[prod(d) N]);
msk = reshape(msk,[prod(d) N]);

msk = sum(msk,2)==N;
nm  = nnz(msk);

if nargout==1
    % Only return number of voxels
    return
end

% Mask images
for n=1:N
    f(~msk,n) = NaN;
end

% Make images have simillar means
for n=1:N
    a      = 512/mean(f(msk,n));
    f(:,n) = a*f(:,n);
end

% Label images using k-means
labels = label_data(f,K,d);

% Generate estimates of MoG parameters
mom        = spm_SuffStats(f,labels);
[wp,mn,vr] = spm_GaussiansFromSuffStats(mom);
wp         = wp/sum(wp); 

% Sort estimates according to mean values of first image
[~,ix] = sort(mn(1,:),2);

mn = mn(:,ix);
vr = vr(:,:,ix);
wp = wp(:,ix);   
%==========================================================================

%==========================================================================
function nlabels = label_data(f,K,d)
% w = warning('query','last')
warning('off','stats:kmeans:MissingDataRemoved')
warning('off','stats:kmeans:FailedToConvergeRep')

opts = statset('MaxIter',500);

labels = kmeans(f,single(K),...
                'Distance','cityblock',...
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