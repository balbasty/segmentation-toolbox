function [mu,a,L] = spm_shoot_blur_wp_load(pth_obj,pth_H,pth_gr,mu,prm,iter)
% A function for blurring ("smoothing") tissue probability maps
% FORMAT [sig,a_new] = spm_shoot_blur(t,prm,its,sig)
%     t   - sufficient statistics
%     prm - regularisation parameters (1,1,1, 0.01,0.02,1)
%     its - max no. iterations (12)
%     sig - optional starting estimates
%
%     sig - "smoothed" average
%     a   - parameters
%
% The core of this procedure is described in:
%     John Ashburner & Karl J. Friston.
%     "Computing Average Shaped Tissue Probability Templates"
%     NeuroImage, In Press, Accepted Manuscript, Available online 24 December 2008
%
% However, there is an additional modification such that the the null space
% of the parameters is rotated out.
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2009)

% John Ashburner
% $Id: spm_shoot_blur.m 5485 2013-05-09 15:51:24Z john $

ll    = 0;
for m=1:numel(pth_obj)
    S = numel(pth_obj{m}); 
    for s=1:S
        obj = load(pth_obj{m}{s},'ll_template'); 
        ll  = ll + obj.ll_template;
    end
end
clear obj

d      = [size(mu),1,1,1];
rits   = [1 1]; % No. cycles and no. relaxation iterations
        
% Only d(4)-1 fields need to be estimated because sum(a,4) = 0.  This matrix
% is used to rotate out the null space
R     = null(ones(1,d(4)));

% Initial starting estimates
a = zeros([d(1:3),d(4)-1],'single');
for z=1:d(3), % Loop over planes
    sz = mu(:,:,z,:);
%     sz = min(max(sig,0),1); % FIX?
    sz = min(max(sz,0),1);
    sz(~isfinite(sz)) = 1/d(4);
    sz = squeeze(log(double(sz*(1-d(4)*1e-3)+1e-3)));
    for j1=1:(d(4)-1),
        az = zeros(d(1:2));
        for j2=1:d(4),
            az = az + R(j2,j1)*sz(:,:,j2); % Note the rotation
        end
        a(:,:,z,j1) = az;
    end
    clear sz az
end

Nii = nifti(pth_gr);
gr  = single(Nii.dat(:,:,:,:));

Nii = nifti(pth_H);
W   = single(Nii.dat(:,:,:,:));
clear Nii
 
% ss1 and ss2 are for examining how close the 1st derivatives are to zero.
% At convergence, the derivatives from the likelihood term should match those
% from the prior (regularisation) term.
% ss1 = sum(sum(sum(sum(gr.^2))));
gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
ll1 = 0.5*sum(sum(sum(sum(gr1.*a)))); % -ve log probability of the prior term
gr  = gr + gr1;                       % Combine the derivatives of the two terms
ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence
mx  = max(max(max(sum(gr.^2,4))));

fprintf('%2d %8.4f %8.4f %8.4f %g\n',iter,ll/prod(d(1:3)),ll1/prod(d(1:3)),(ll+ll1)/prod(d(1:3)),(ss2)/prod(d(1:3)));

reg = double(0.01*sqrt(mx)*d(4));
%reg = double(0.1*sqrt(ss2/prod(d(1:3))));
a   = a - spm_field(W,gr,[prm(1:3) prm(4)+reg prm(5:6) rits]); % Gauss-Newton update  
    
mu = sftmax(a,R);
L  = -(ll + ll1);
%________________________________________________________

%________________________________________________________
function sig = sftmax(a,R,log_wp)
% Softmax function
if nargin<3, log_wp = 0; end

d     = [size(a) 1 1 1];
sig   = zeros([d(1:3),d(4)+1],'single');
trunc = log(realmax('single')*(1-eps('single'))/(d(4)+1));

for j=1:size(a,3), % Loop over planes

    % Rotate the null-space back in to the data
    aj  = double(reshape(a(:,:,j,:),[d(1:2),d(4)]));
    sj  = zeros([d(1:2),d(4)+1]);
    for j1=1:d(4)+1,
        sj(:,:,j1) = 0;
        for j2=1:d(4),
            sj(:,:,j1) = sj(:,:,j1) + R(j1,j2)*aj(:,:,j2);
        end
    end

    % Compute safe softmax
    sj           = bsxfun(@plus,sj,log_wp);
    mx_sj        = max(sj,[],3);
    sj           = exp(bsxfun(@minus,sj,mx_sj));
    s            = sum(sj,3);
    sig(:,:,j,:) = single(bsxfun(@rdivide,sj,s));
end
%________________________________________________________