function V = read_images(im,pars)
pth      = im{1};
S        = im{2};
modality = im{3};

do_rem_corrupted = pars.preproc.do_rem_corrupted;
tol_dist         = pars.preproc.tol_dist;
tol_vx           = pars.preproc.tol_vx;
verbose          = pars.preproc.verbose;
    
if do_rem_corrupted, num_workers = Inf;
else                 num_workers = 0;
end

folder    = dir(pth);  % folder with subfolders containing multi-channel data of subjects
folder    = folder(3:end);
dirflag   = [folder.isdir];
subfolder = folder(dirflag); % subfolders (S1,S2,...)
S1        = numel(subfolder);
if S>S1
    S = S1;
end
files     = dir(fullfile(pth,subfolder(1).name,'*.nii'));
N         = numel(files);

V    = cell(1,S);
sd   = zeros(N,S);
v    = zeros(N,S);
sint = zeros(N,S);
vx   = zeros(N,S);

spm_parfor('manage_parpool',num_workers);
parfor (s=1:S,num_workers) 
    folder = fullfile(pth,subfolder(s).name);
    files  = dir(fullfile(folder,'*.nii'));
    for n=1:N
        V{s}(n) = spm_vol(fullfile(folder,files(n).name));
        
        if do_rem_corrupted            
            [sd(n,s),v(n,s),sint(n,s),vx(n,s)] = compute_img_stats(V{s}(n).fname,modality);                     
        end
    end        
end 

if do_rem_corrupted
    % Standardise the data (zero mean and unit variance)
    X = [sd(:) v(:) sint(:)];
    X = bsxfun(@minus,X,mean(X));
    X = bsxfun(@rdivide,X,sqrt(var(X)));            
    
    ix = zeros(1,S);
    f  = cell(1,S);
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
                if DM>tol_dist || vx(n,s)>tol_vx                               
                    ix(s) = s;
                    f{s}  = V{s}(n).fname;
                end
                cnt = cnt + 1;
            end
        end
    end

    % Remove 'outliers'
    ix(ix==0) = [];
    V(ix)     = [];
    S         = S - numel(ix);

    if numel(ix)>0
        f = f(~cellfun('isempty',f));

        if verbose, spm_check_registration(char(f')); end

        for s=1:numel(f)
           disp(['Removing image ' f{s}]) ;
        end
        fprintf('%d subjects remaining\n',S)
    end    
end

fprintf('Loaded data from %d subject(s) having %d channel(s) each\n',S,N); 
%==========================================================================

%==========================================================================
function [sd,v,sint,vx] = compute_img_stats(pth,modality)

% Get image and mask
Nii     = nifti(pth);
f       = single(Nii.dat(:,:,:));
msk     = msk_modality(f,modality);
f(~msk) = NaN;
clear msk

% Calculate sum of image intensities
sint = nansum(f(:));

% Calculate image gradients
vx              = vxsize(Nii.mat);
[gx,gy,gz]      = spm_imcalc('grad',f,vx);
f(~isfinite(f)) = [];

vx = round(vx);
vx = max(vx,1);
vx = max(vx)/min(vx);

% Estimate FWHM
fwhm = sqrt(4*log(2))*sum(abs(f))./[sum(abs(gx(~isfinite(gx)))) sum(abs(gy(~isfinite(gy)))) sum(abs(gz(~isfinite(gz))))];
fwhm(~isfinite(fwhm)) = 0;

% Estimate noise standard deviation
sc = prod([spm_smoothkern(fwhm(1),0) spm_smoothkern(fwhm(2),0) spm_smoothkern(fwhm(3),0)]);
sc = sc/2;
sc = min(sc,1);
sd = sqrt(sum(f.^2)/(numel(f)*sc));
clear f

% Calculate gradient magnitude
gm                = sqrt(gx.^2 + gy.^2 + gz.^2);
gm(~isfinite(gm)) = [];
clear gx gy gz

% Fit histogram to gradient magnitude (GM)
[h,x] = hist(gm,10000);
p     = h/numel(gm)/mean(diff(x));
v     = sum(p.*x.^2)/sum(p);
%==========================================================================

%==========================================================================
function val = dist_Mahalanobis(x,mu,C)
val = sqrt((x - mu)'*(C\(x - mu))); % Mahalanobis distance
%==========================================================================