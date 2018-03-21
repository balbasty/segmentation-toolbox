function [mu,L] = spm_shoot_blur_wp(obj,mu,prm,iter,its,verbose)
if nargin<6, verbose = false; end

rits = [1 1]; % No. cycles and no. relaxation iterations

M = numel(obj);
d1 = [size(mu),1,1,1];
        
% Only d(4)-1 fields need to be estimated because sum(a,4) = 0.  This matrix
% is used to rotate out the null space
R = null(ones(1,d1(4)));

a = zeros([d1(1:3),d1(4)-1],'single');
for z=1:d1(3), % Loop over planes
    sz = mu(:,:,z,:);
    for j1=1:(d1(4)-1),
        az = zeros(d1(1:2));
        for j2=1:d1(4),
            az = az + R(j2,j1)*sz(:,:,j2); % Note the rotation
        end
        a(:,:,z,j1) = az;
    end
    clear sz az
end

for i=1:its,
    W  = zeros([d1(1:3) round(((d1(4)-1)*d1(4))/2)],'single'); % 2nd derivatives
    gr = zeros([d1(1:3),d1(4)-1],'single');                   % 1st derivatives    
    ll = 0;
    
    tot_subj = 0;
    for m=1:M        
        S        = numel(obj{m});    
        tot_subj = tot_subj + S;
                
        for s1=1:S % Loop over subjects

            [t,s,lwp,rng0,d0] = load_from_obj(obj,m,s1);                        
            lwp1              = reshape(lwp,1,1,1,d1(4));
            
            % Compute gradients and Hessian
            lls = 0;            
            for z=1:d0(3), % Loop over planes

                a1 = a(rng0{1},rng0{2},rng0{3}(z),:);       

                % Compute softmax for this plane
                mu = double(reshape(sftmax(a1,R,lwp),[d0(1:2),d1(4)]));

                % -ve log likelihood of the likelihood
%                 ll0 = sum(sum(sum(log(mu).*reshape(t(:,:,z,:),[d0(1:2),d1(4)]),3).*s(:,:,z)));
                                
                % log-likelihood (using log-sum-exp)
                sm0 = bsxfun(@plus,rotate_back(a1,R),lwp1);
                sm1 = sum(t(:,:,z,:).*sm0,4);
                sm2 = logsumexp(sm0,4).*sum(t(:,:,z,:),4);
                ll0 = sum(sum(sm1 - sm2));                                
                
                lls = lls - ll0;
                
                if ~isfinite(lls)
                    warning('~isfinite(lls)');
                end
                
                % Compute first derivatives (d(4)-1) x 1 
                grz = bsxfun(@times,mu,s(:,:,z)) - double(reshape(t(:,:,z,:),[d0(1:2),d1(4)]));
                for j1=1:(d1(4)-1),
                    gr1 = zeros([d0(1:2) 1 d1(4) - 1],'single');
                    for j2=1:d1(4),
                        gr1(:,:,1,j1) = gr1(:,:,1,j1) + R(j2,j1)*grz(:,:,j2); % Note the rotation
                    end
                    gr(rng0{1},rng0{2},rng0{3}(z),j1) = gr(rng0{1},rng0{2},rng0{3}(z),j1) + gr1(:,:,1,j1);
                end

                % Compute d(4) x d(4) matrix of second derivatives at each voxel.
                % These should be positive definate, but rounding errors may prevent this.
                % Regularisation is included to enforce +ve definateness.
                wz = zeros([d0(1:2),d1(4),d1(4)]);
                for j1=1:d1(4),
                    wz(:,:,j1,j1) =   (1-mu(:,:,j1)).*mu(:,:,j1).*s(:,:,z);
                    for j2=1:(j1-1),
                        wz(:,:,j1,j2) = -mu(:,:,j1) .*mu(:,:,j2).*s(:,:,z);
                        wz(:,:,j2,j1) = wz(:,:,j1,j2);
                    end
                end

                % First step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
                % by R'*W*R
                wz1 = zeros([d0(1:2),d1(4),d1(4)-1]);
                for j1=1:d1(4),
                    for j2=1:(d1(4)-1),
                        tmp = zeros(d0(1:2));
                        for j3=1:d1(4),
                            tmp = tmp + wz(:,:,j1,j3)*R(j3,j2);
                        end
                        wz1(:,:,j1,j2) = tmp;
                    end
                end

                % Second step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
                % by R'*W*R
                wz = zeros([d0(1:2),d1(4)-1,d1(4)-1]);
                for j1=1:(d1(4)-1),
                    for j2=1:(d1(4)-1),
                        tmp = zeros(d0(1:2));
                        for j3=1:d1(4),
                            tmp = tmp + R(j3,j1)*wz1(:,:,j3,j2);
                        end
                        wz(:,:,j1,j2) = tmp;
                    end
                end

                % First pull out the diagonal of the 2nd derivs
                for j1=1:d1(4)-1,
                    W(rng0{1},rng0{2},rng0{3}(z),j1) = W(rng0{1},rng0{2},rng0{3}(z),j1) + wz(:,:,j1,j1) + 1e-4*d(4)^2;
                end

                % Then pull out the off diagonal parts (note that matrices are symmetric)
                jj = d1(4);
                for j1=1:d1(4)-1,
                   for j2=(j1+1):(d1(4)-1),
                       W(rng0{1},rng0{2},rng0{3}(z),jj) = W(rng0{1},rng0{2},rng0{3}(z),jj) + wz(:,:,j2,j1);
                       jj = jj+1;
                   end
                end
            end
            ll = ll + lls;
        end        
    end
    
    % ss1 and ss2 are for examining how close the 1st derivatives are to zero.
    % At convergence, the derivatives from the likelihood term should match those
    % from the prior (regularisation) term.
    ss1 = sum(sum(sum(sum(gr.^2))));    
    gr1 = spm_field('vel2mom',a,prm);     % 1st derivative of the prior term
    ll1 = 0.5*sum(sum(sum(sum(gr1.*a)))); % -ve log probability of the prior term
    gr  = gr + gr1;                       % Combine the derivatives of the two terms
    ss2 = sum(sum(sum(sum(gr.^2))));      % This should approach zero at convergence
    
    if verbose
        fprintf('%2d | %8.7f %8.7f %8.7f %8.7f\n',iter,ll/prod(d1(1:3)),ll1/prod(d1(1:3)),(ll + ll1)/prod(d1(1:3)),(ss2)/prod(d1(1:3)));
    end
    
    a = a - spm_field(W,gr,[prm(1:3) prm(4) prm(5:6) rits]); % Gauss-Newton update  
    
    if ss2/ss1<1e-4, break; end        % Converged?
end
mu = rotate_back(a,R);
L  = -(ll + ll1);
%==========================================================================

%==========================================================================
function sig = sftmax(a,R,lwp)
% Softmax function
if nargin<3, lwp = 0; end

d   = [size(a) 1 1 1];
sig = zeros([d(1:3),d(4)+1],'single');

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
    sj           = bsxfun(@plus,sj,lwp);
    mx_sj        = max(sj,[],3);
    sj           = exp(bsxfun(@minus,sj,mx_sj));
    s            = sum(sj,3);
    sig(:,:,j,:) = single(bsxfun(@rdivide,sj,s));
end
%==========================================================================

%==========================================================================
function [t,s,lwp,rng0,d1] = load_from_obj(obj,m,s1)
K   = numel(obj{m}{s1}.pth_resp);
Nii = nifti(obj{m}{s1}.pth_resp{1});
dm  = size(Nii.dat(:,:,:));
if numel(dm)==2, dm(3) = 1; end
dm  = [dm K];

t          = zeros(dm,'single');
t(:,:,:,1) = Nii.dat(:,:,:);
for k=2:K
    Nii        = nifti(obj{m}{s1}.pth_resp{k});
    t(:,:,:,k) = Nii.dat(:,:,:);
end
t   = max(t,eps('single')*1000);
s   = sum(t,4);
lwp = reshape(log(obj{m}{s1}.wp),1,1,dm(4));

bb_push  = obj{m}{s1}.bb_push;
rng0     = cell(1,3);
rng0{1}  = bb_push(1,1):bb_push(1,2);
rng0{2}  = bb_push(2,1):bb_push(2,2);
rng0{3}  = bb_push(3,1):bb_push(3,2);            

d1 = dm(1:3);
%==========================================================================

%==========================================================================
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
%==========================================================================
