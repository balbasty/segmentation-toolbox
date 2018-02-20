function [mu,L] = spm_shoot_blur_wp_load(obj,mu,prm,iter,verbose)
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

if nargin<5, verbose = false; end

d    = [size(mu),1,1,1];
rits = [1 1]; % No. cycles and no. relaxation iterations

tot_subj = 0;
ll       = 0;
W        = zeros([d(1:3) round(((d(4)-1)*d(4))/2)],'single'); % 2nd derivatives
gr       = zeros([d(1:3),d(4)-1],'single');                   % 1st derivatives
for m=1:numel(obj)
    S = numel(obj{m}); 
    for s=1:S
        if obj{m}{s}.status==0
            bb_push = obj{m}{s}.bb_push;
            rngx    = bb_push(1,1):bb_push(1,2);
            rngy    = bb_push(2,1):bb_push(2,2);
            rngz    = bb_push(3,1):bb_push(3,2);
            
            for k=1:size(gr,4)
                Nii                  = nifti(obj{m}{s}.pth_gr{k});
                gr(rngx,rngy,rngz,k) = gr(rngx,rngy,rngz,k) + single(Nii.dat(:,:,:));
            end
            
            for k=1:size(W,4)
                Nii                 = nifti(obj{m}{s}.pth_H{k});
                W(rngx,rngy,rngz,k) = W(rngx,rngy,rngz,k) + single(Nii.dat(:,:,:));
            end
            
            ll       = ll + obj{m}{s}.ll_template;            
            tot_subj = tot_subj + 1;            
        end
    end
end
clear obj Nii
        
prm(4) = prm(4) + tot_subj*d(4)*1e-6;

% Only d(4)-1 fields need to be estimated because sum(a,4) = 0.  This matrix
% is used to rotate out the null space
R = null(ones(1,d(4)));

% Initial starting estimates
a = zeros([d(1:3),d(4)-1],'single');
for z=1:d(3), % Loop over planes
    sz = mu(:,:,z,:);
    for j1=1:(d(4)-1),
        az = zeros(d(1:2));
        for j2=1:d(4),
            az = az + R(j2,j1)*sz(:,:,j2); % Note the rotation
        end
        a(:,:,z,j1) = az;
    end
    clear sz az
end

% ss1 and ss2 are for examining how close the 1st derivatives are to zero.
% At convergence, the derivatives from the likelihood term should match those
% from the prior (regularisation) term.
% ss1 = sum(sum(sum(sum(gr.^2))));
gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
ll1 = 0.5*sum(sum(sum(sum(gr1.*a)))); % -ve log probability of the prior term
gr  = gr + gr1;                       % Combine the derivatives of the two terms
ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence
mx  = max(max(max(sum(gr.^2,4))));

if verbose
    fprintf('%2d %8.7f %8.7f %8.7f %8.7f\n',iter,ll/prod(d(1:3)),ll1/prod(d(1:3)),(ll + ll1)/prod(d(1:3)),(ss2)/prod(d(1:3)));
end

reg = double(0.01*sqrt(mx)*d(4));
%reg = double(0.1*sqrt(ss2/prod(d(1:3))));
a   = a - spm_field(W,gr,[prm(1:3) prm(4)+reg prm(5:6) rits]); % Gauss-Newton update  
    
mu = rotate_back(a,R);
L  = -(ll + ll1);
%________________________________________________________

%________________________________________________________
function sig = rotate_back(a,R)
d   = [size(a) 1 1 1];
sig = zeros([d(1:3),d(4) + 1],'single');

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

    sig(:,:,j,:) = sj;
end
%________________________________________________________