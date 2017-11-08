function [V,S] = rem_corrupted_ims(V,verbose)
if nargin<2, verbose = 0; end

% Calculate some image statistics that will be used as indicators if images
% are corrupted or not
S  = numel(V);
N  = numel(V{1});
sc = zeros(1,S*N);
v  = zeros(1,S*N);
for s=1:S
    for n=1:N
        [~,sc(s),~,v(s)] = compute_img_stats(V{s}(n).fname);
    end
end

% Make data have zero mean and unit variance
X = [sc', v'];
X = bsxfun(@minus,X,mean(X,1));
X = bsxfun(@rdivide,X,var(X));
clear sc v

% Fit a Gaussian to the image statistics
mu = mean(X)';
C  = cov(X);

% dm  = zeros(1,S);
% for s=1:S
%     x     = [X(s,1),X(s,2)]';
%     dm(s) = sqrt((x - mu)'*(C\(x - mu))); % Mahalanobis distance
% end
% scatter3(X(:,1)',X(:,2)',dm)

% Remove images based on Mahalanobis distance
tol = 4; % Seems, empirically, to be a descent value
ix  = [];
f   = {};
cnt = 1;
for s=1:S
    for n=1:N
        x  = [X(cnt,1),X(cnt,2)]';
        DM = sqrt((x - mu)'*(C\(x - mu))); % Mahalanobis distance
        if DM>tol
            disp(['Removed CT image: ' V{s}(n).fname])
            f{end + 1} = V{s}(n).fname;
            ix         = [ix,s];
        end
        cnt = cnt +1;
    end
end
V(ix) = [];
S     = S - numel(ix);

if numel(ix)>0
    if verbose, spm_check_registration(char(f)); end

    fprintf('%d subjects remaining\n',S)
end

% % Fit Gaussians to each image statistic
% msc   = mean(sc);
% sigsc = sqrt(var(sc));
% 
% mv   = mean(v);
% sigv = sqrt(var(v));
% 
% % Define thresholds based on the distance from the mean (measured in SDs)
% sds   = 3; % This seems empirically to be a good value..
% tolsc = msc + sds*sigsc;
% tolv  = mv  + sds*sigv;
% 
% % Remove images based on the above thresholds
% ix   = [];
% f    = {};
% for s=1:S
%     if sc(s)>tolsc && v(s)>tolv
%         disp(['Removed CT image: ' V{s}.fname])
%         f{end + 1} = V{s}.fname;
%         ix         = [ix,s];
%     end
% end
% V(ix) = [];
% S     = S - numel(ix);
%==========================================================================

%==========================================================================
function [fwhm,sc,sd,v] = compute_img_stats(pth)

% Calculate image gradients
Nii        = nifti(pth);
f          = Nii.dat(:,:,:);
msk        = get_msk(f);
f(~msk)    = NaN;
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