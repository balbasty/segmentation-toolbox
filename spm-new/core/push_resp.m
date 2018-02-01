function ll = push_resp(obj,tpm,bb,vx)
if nargin<3, bb = NaN(2,3);   end % Default to TPM bounding box
if nargin<4, vx = NaN;        end % Default to TPM voxel size

lkp          = obj.lkp;
wp           = obj.wp;
pth_template = obj.pth_template;
Kb           = max(lkp);
N            = numel(obj.image);
modality     = obj.modality;

% Read essentials from tpm (it will be cleared later)
%--------------------------------------------------------------------------
M1 = tpm.M;

% Define orientation and field of view of any "normalised" space
% data that may be generated (wc*.nii, mwc*.nii, rc*.nii & y_*.nii).
%--------------------------------------------------------------------------
if any(isfinite(bb(:))) || any(isfinite(vx))
    % If a bounding box is supplied, combine this with the closest
    % bounding box derived from the dimensions and orientations of
    % the tissue priors.
    [bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
    bb(~isfinite(bb)) = bb1(~isfinite(bb));
    if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end
    bb(1,:) = vx*round(bb(1,:)/vx);
    bb(2,:) = vx*round(bb(2,:)/vx);
    odim    = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;

    mm  = [[bb(1,1) bb(1,2) bb(1,3)
            bb(2,1) bb(1,2) bb(1,3)
            bb(1,1) bb(2,2) bb(1,3)
            bb(2,1) bb(2,2) bb(1,3)
            bb(1,1) bb(1,2) bb(2,3)
            bb(2,1) bb(1,2) bb(2,3)
            bb(1,1) bb(2,2) bb(2,3)
            bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];
    vx3 = [[1       1       1
            odim(1) 1       1
            1       odim(2) 1
            odim(1) odim(2) 1
            1       1       odim(3)
            odim(1) 1       odim(3)
            1       odim(2) odim(3)
            odim(1) odim(2) odim(3)]'; ones(1,8)];
    mat    = mm/vx3;
else
    % Use the actual dimensions and orientations of
    % the tissue priors.
    odim = tpm.V(1).dim;
    mat  = tpm.V(1).mat;
end

%--------------------------------------------------------------------------
d         = obj.image(1).dim(1:3);
[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3        = 1:d(3);

%--------------------------------------------------------------------------
chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N
    d3         = [size(obj.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = obj.Tbias{n};

    % Need to fix writing of bias fields or bias corrected images, when the data used are 4D.    
    chan(n).ind = obj.image(n).n;
end

%--------------------------------------------------------------------------
prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Nii     = nifti(obj.pth_vel);
Twarp   = single(Nii.dat(:,:,:,:)); 
Coef{1} = spm_bsplinc(Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(Twarp(:,:,:,3),prm);
M       = M1\obj.Affine*obj.image(1).mat;
clear Nii Twarp

%--------------------------------------------------------------------------
y = zeros([obj.image(1).dim(1:3),3],'single');
Q = zeros([d(1:3),Kb],'single');
for z=1:length(x3)
    % Bias corrected image
    fN  = cell(1,N);
    bfN = cell(1,N);
    for n=1:N
        fN{n} = spm_sample_vol(obj.image(n),x1,x2,o*x3(z),0);
        fN{n} = fN{n}(:);
        
        if strcmp(modality,'MRI')  
            msk         = isfinite(fN{n}) & (fN{n}~=0);
            fN{n}(~msk) = NaN;            
        elseif strcmp(modality,'CT')
            msk         = isfinite(fN{n}) & (fN{n}~=min(fN{n})) & (fN{n}~=0) & (fN{n}~=-3024) & (fN{n}~=-1500);
            fN{n}(~msk) = NaN;
            fN{n}       = fN{n} + 1000;
        end
      
        bfN{n} = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bfN{n} = bfN{n}(:);      
    end
    clear f 
    
    % Compute the deformation (mapping voxels in image to voxels in TPM)
    [t1,t2,t3] = make_deformation(Coef,z,obj.MT,prm,x1,x2,x3,M);

    y(:,:,z,1) = t1;
    y(:,:,z,2) = t2;
    y(:,:,z,3) = t3;    
    
    % Parametric representation of intensity distributions                
    b = spm_sample_priors8(tpm,t1,t2,t3);
    clear t1 t2 t3
    
    B = zeros([prod(d(1:2)) Kb]);
    for k1 = 1:Kb, B(:,k1) = b{k1}(:); end
    clear b
    
    q  = zeros([d(1:2) Kb]);  
    q1 = latent(fN,bfN,obj.mg,obj.gmm,B,lkp,obj.wp,[]);
    clear B fN bfN
    
    q1 = reshape(q1,[d(1:2),numel(obj.mg)]);
    for k1=1:Kb
        tmp       = sum(q1(:,:,lkp==k1),3);
        q(:,:,k1) = tmp;
    end 
    clear q1
    
    Q(:,:,z,:) = single(reshape(q,[d(1:2),1,Kb]));
    clear q
end
clear tpm x1 x2 x3 o Coeff chan

% Adjust stuff so that warped data (and deformations) have the
% desired bounding box and voxel sizes, instead of being the same
% as those of the tissue probability maps.
%--------------------------------------------------------------------------
M = mat\M1;
for i=1:size(y,3)
    t1         = y(:,:,i,1);
    t2         = y(:,:,i,2);
    t3         = y(:,:,i,3);
    y(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
    y(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
    y(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
end
M1 = mat;
d1 = odim;
clear t1 t2 t3 

% Push responsibilities to template space
%--------------------------------------------------------------------------
t = zeros([d1 Kb],'single');
for k1=1:Kb
    t(:,:,:,k1) = spm_diffeo('push',Q(:,:,:,k1),y,d1(1:3));

    fname  = obj.pth_resp{k1};
    if (exist(fname,'file')==2), delete(fname); end
    
    Nii      = nifti;
    Nii.dat  = file_array(fname,...
                        d1,...
                        [spm_type('float32') spm_platform('bigend')],...
                        0,1,0);
    Nii.mat  = M1;
    Nii.mat0 = M1;
    Nii.descrip = ['Pushed responsibilities ' num2str(k1)];
    create(Nii);
    
    Nii.dat(:,:,:) = t(:,:,:,k1);    
end
clear Q y Nii

% Compute template gradient and Hessian
%--------------------------------------------------------------------------
Nii = nifti(pth_template);
mu  = single(Nii.dat(:,:,:,:));
d   = [size(mu),1,1,1];
clear Nii

% Only d(4)-1 fields need to be estimated because sum(a,4) = 0.  This matrix
% is used to rotate out the null space
R = null(ones(1,d(4)));

% Initial starting estimates
a = zeros([d(1:3),d(4)-1],'single');
for z=1:d(3), % Loop over planes
    sz = mu(:,:,z,:);
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

W   = zeros([d(1:3) round(((d(4)-1)*d(4))/2)],'single'); % 2nd derivatives
gr  = zeros([d(1:3),d(4)-1],'single');                   % 1st derivatives
dgr = size(gr);

% Read a responsibility and log(wp)      
lwp = reshape(log(wp),1,1,d(4));

% Re-organise sufficient statistics to a form that is easier to work with
t = max(t,eps('single')*1000);
s = sum(t,4);
for k=1:d(4)
    t(:,:,:,k) = t(:,:,:,k)./s;
end

% Compute gradients and Hessian
ll = 0;
for z=1:d(3), % Loop over planes

    % Compute softmax for this plane
    mu = double(reshape(sftmax(a(:,:,z,:),R,lwp),[d(1:2),d(4)]));

    % -ve log likelihood of the likelihood
    ll  = ll - sum(sum(sum(log(mu).*reshape(t(:,:,z,:),[d(1:2),d(4)]),3).*s(:,:,z)));

    % Compute first derivatives (d(4)-1) x 1 
    grz = mu - double(reshape(t(:,:,z,:),[d(1:2),d(4)]));
    for j1=1:(d(4)-1),
        gr1 = zeros([dgr(1:2) 1 dgr(4)],'single');
        for j2=1:d(4),
            gr1(:,:,1,j1) = gr1(:,:,1,j1) + R(j2,j1)*grz(:,:,j2); % Note the rotation
        end
        gr(:,:,z,j1) = gr(:,:,z,j1) + gr1(:,:,1,j1).*s(:,:,z);
    end

    % Compute d(4) x d(4) matrix of second derivatives at each voxel.
    % These should be positive definate, but rounding errors may prevent this.
    % Regularisation is included to enforce +ve definateness.
    wz = zeros([d(1:2),d(4),d(4)]);
    for j1=1:d(4),
        wz(:,:,j1,j1) =   (1-mu(:,:,j1)).*mu(:,:,j1).*s(:,:,z);
        for j2=1:(j1-1),
            wz(:,:,j1,j2) = -mu(:,:,j1) .*mu(:,:,j2).*s(:,:,z);
            wz(:,:,j2,j1) = wz(:,:,j1,j2);
        end
    end

    % First step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
    % by R'*W*R
    wz1 = zeros([d(1:2),d(4),d(4)-1]);
    for j1=1:d(4),
        for j2=1:(d(4)-1),
            tmp = zeros(d(1:2));
            for j3=1:d(4),
                tmp = tmp + wz(:,:,j1,j3)*R(j3,j2);
            end
            wz1(:,:,j1,j2) = tmp;
        end
    end

    % Second step of rotating 2nd derivatives to (d(4)-1) x (d(4)-1)
    % by R'*W*R
    wz = zeros([d(1:2),d(4)-1,d(4)-1]);
    for j1=1:(d(4)-1),
        for j2=1:(d(4)-1),
            tmp = zeros(d(1:2));
            for j3=1:d(4),
                tmp = tmp + R(j3,j1)*wz1(:,:,j3,j2);
            end
            wz(:,:,j1,j2) = tmp;
        end
    end

    % First pull out the diagonal of the 2nd derivs
    for j1=1:d(4)-1,
        W(:,:,z,j1) = W(:,:,z,j1) + wz(:,:,j1,j1);% + maxs*sqrt(eps('single'))*d(4)^2;
    end

    % Then pull out the off diagonal parts (note that matrices are symmetric)
    jj = d(4);
    for j1=1:d(4)-1,
       for j2=(j1+1):(d(4)-1),
           W(:,:,z,jj) = W(:,:,z,jj) + wz(:,:,j2,j1);
           jj = jj+1;
       end
    end
end

% Save gradients and Hessians
Nii              = nifti(obj.pth_gr);
Nii.dat(:,:,:,:) = Nii.dat(:,:,:,:) + gr;

Nii              = nifti(obj.pth_H);
Nii.dat(:,:,:,:) = Nii.dat(:,:,:,:) + W;
%==========================================================================

%==========================================================================
function [x1,y1,z1] = make_deformation(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%==========================================================================

%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end
return;
%==========================================================================

%==========================================================================
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
%==========================================================================