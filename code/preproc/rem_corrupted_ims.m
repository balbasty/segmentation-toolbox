function [V,labels] = rem_corrupted_ims(V,labels,num_workers,descrip,verbose)
if nargin<5, verbose = 0; end

% Calculate some image statistics that will be used as indicators if images
% are corrupted or not
S    = numel(V);
N    = numel(V{1});
sc   = zeros(1,S*N);
sd   = zeros(1,S*N);
v    = zeros(1,S*N);
sint = zeros(1,S*N);
parfor (s=1:S,num_workers)
    for n=1:N
        [~,sc(s),sd(s),v(s),sint(s)] = compute_img_stats(V{s}(n).fname,descrip);
    end
end

% Standardise the data (zero mean and unit variance)
X = [sc',sd',v',sint'];
X = bsxfun(@minus,X,mean(X));
X = bsxfun(@rdivide,X,sqrt(var(X)));
clear sc sd v sint

tol = 4; % Seems, empirically, to be a descent value
ix  = zeros(1,S);
f   = cell(1,S);
for i=1:size(X,2)

    % Fit a Gaussian to the data
    mu = mean(X(:,i))';
    C  = cov(X(:,i));

    if verbose==2
        dm  = zeros(1,size(X,1));
        for s=1:size(X,1)
            x     = X(s,i)';
            dm(s) = dist_Mahalanobis(x,mu,C);
        end
        figure;
        scatter(X(:,i)',dm)
        drawnow
    end
        
    cnt = 1;
    for s=1:S
        for n=1:N
            x  = X(cnt,i)';
            DM = dist_Mahalanobis(x,mu,C);
            if DM>tol                                
                ix(s) = s;
                f{s}  = V{s}(n).fname;
            end
            cnt = cnt + 1;
        end
    end
end

% Remove 'outliers'
ix(ix==0)  = [];
V(ix)      = [];
labels(ix) = [];
S          = S - numel(ix);

if numel(ix)>0
    f = f(~cellfun('isempty',f));
    
    if verbose, spm_check_registration(char(f')); end

    for s=1:numel(f)
       disp(['Removing image ' f{s}]) ;
    end
    fprintf('%d subjects remaining\n',S)
end
%==========================================================================

%==========================================================================
function val = dist_Mahalanobis(x,mu,C)
val = sqrt((x - mu)'*(C\(x - mu))); % Mahalanobis distance
%==========================================================================

%==========================================================================
function [fwhm,sc,sd,v,sint] = compute_img_stats(pth,descrip)

% Get image and mask
Nii        = nifti(pth);
f          = Nii.dat(:,:,:);
msk        = get_msk(f,descrip);
f(~msk)    = NaN;

% Calculate sum of image intensities
sint = nansum(abs(f(:)));

% Calculate image gradients
vx         = sqrt(sum(Nii.mat(1:3,1:3).^2));
[gx,gy,gz] = grad(f,vx);

% Estimate FWHM
fwhm = sqrt(4*log(2))*sum(abs(f(isfinite(f))))./[sum(abs(gx(isfinite(gx)))) sum(abs(gy(isfinite(gy)))) sum(abs(gz(isfinite(gz))))];
fwhm(~isfinite(fwhm)) = 0;

% Estimate noise standard deviation
sc = prod([spm_smoothkern(fwhm(1),0) spm_smoothkern(fwhm(2),0) spm_smoothkern(fwhm(3),0)]);
sc = sc/2;
sc = min(sc,1);
sd = sqrt(sum(f(isfinite(f)).^2)/(numel(f(isfinite(f)))*sc));

% Fit histogram to gradient magnitude (GM)
gm      = sqrt(gx.^2 + gy.^2 + gz.^2);
[h,x]   = hist(gm(isfinite(gm)),10000);

% Smooth histogram a little bit
fwhmh = 3;
lim   = ceil(2*fwhmh);
krn1  = spm_smoothkern(fwhmh(1),-lim(1):lim(1));
krn1  = krn1/sum(krn1); 
h     = conv2(h,krn1,'same');

% Compute variance from GM histogram
p = h/numel(gm(isfinite(gm)))/mean(diff(x));
v = sum(p.*x.^2)/sum(p);
%==========================================================================
