function [pthx,pthmsk,mat,d] = warp_im_same_size(Px,samp,tempdir,ord)
if nargin < 4, ord = [1 1 1 0 0 0]; end

% Create folders, if doens not exist---------------------------------------
fw = fullfile(tempdir,'warped');
if (exist(fw,'dir') == 0)
    mkdir(fw);
end

fmsk = fullfile(tempdir,'msk');
if (exist(fmsk,'dir') == 0)
    mkdir(fmsk);
end

% Compute average orientation matrix---------------------------------------
N   = numel(Px);
mat = [];
d   = [];
for n=1:N
    Nii = nifti(Px{n});
    mat = cat(3,mat,Nii.mat);
    d   = cat(1,d,size(Nii.dat));
end

[mat,d] = compute_avg_mat(mat,d);

vx       = sqrt(sum(mat(1:3,1:3).^2));
st       = samp./vx;
F        = diag(st);
F(1:4,4) = 1;
mat      = mat*F;
d        = floor(d./st);

% Warp images and create masks---------------------------------------------
pthx   = cell(N,1);
pthmsk = cell(N,1);
for n=1:N 
    Nii  = nifti(Px{n});
    matn = Nii.mat;
    
    % Warp image-----------------------------------------------------------
    x  = single(Nii.dat(:,:,:));            
    x  = warp(x,matn\mat,[],ord,d,[1 1 1]);
    
    % Create mask----------------------------------------------------------
    msk = x~=0 & x~=max(x(:)) & x~=min(x(:)) & isfinite(x) & x~=-3024 & x~=-1500;
    
    % Make images have simillar means--------------------------------------                       
    x = x*(1024/mean(x(msk)));
    
    % Write warped image and mask------------------------------------------
    [~,nam,ext] = fileparts(Nii.dat.fname);
    
    pthx{n} = fullfile(fw,['w' nam ext]);    
    create_nii(pthx{n},x,mat,'float32','X');
    
    pthmsk{n} = fullfile(fmsk,['msk' nam ext]);
    create_nii(pthmsk{n},msk,mat,'uint8-le','Mask');
end
%==========================================================================

%==========================================================================
function [M_avg,d] = compute_avg_mat(Mat0,dims)
% Compute an average voxel-to-world mapping and suitable dimensions
% FORMAT [M_avg,d] = compute_avg_mat(Mat0,dims)
% Mat0  - array of matrices (4x4xN)
% dims  - image dimensions (Nx3)
% M_avg - voxel-to-world mapping
% d     - dimensions for average image
%

% Rigid-body matrices computed from exp(p(1)*B(:,:,1)+p(2)+B(:,:,2)...)
%-----------------------------------------------------------------------
B = se3_basis;

% Find combination of 90 degree rotations and flips that brings all
% the matrices closest to axial
%-----------------------------------------------------------------------
Matrices = Mat0;
pmatrix  = [1,2,3; 2,1,3; 3,1,2; 3,2,1; 1,3,2; 2,3,1];
for i=1:size(Matrices,3)
    vx    = sqrt(sum(Matrices(1:3,1:3,i).^2));
    tmp   = Matrices(:,:,i)/diag([vx 1]);
    R     = tmp(1:3,1:3);
    minss = Inf;
    minR  = eye(3);
    for i1=1:6,
        R1 = zeros(3);
        R1(pmatrix(i1,1),1)=1;
        R1(pmatrix(i1,2),2)=1;
        R1(pmatrix(i1,3),3)=1;
        for i2=0:7,
            F  = diag([bitand(i2,1)*2-1, bitand(i2,2)-1, bitand(i2,4)/2-1]);
            R2 = F*R1;
            ss = sum(sum((R/R2-eye(3)).^2));
            if ss<minss,
                minss = ss;
                minR  = R2;
            end
        end
    end
    rdim = abs(minR*dims(i,:)');
    R2   = inv(minR);
    minR = [R2 R2*((sum(R2,1)'-1)/2.*(rdim+1)); 0 0 0 1];
    Matrices(:,:,i) = Matrices(:,:,i)*minR;
end

% Average of these matrices
%-----------------------------------------------------------------------
M_avg = spm_meanm(Matrices);

% If average involves shears, then find the closest matrix that does not
% require them
%-----------------------------------------------------------------------
p = spm_imatrix(M_avg);
if sum(p(10:12).^2)>1e-8,

    % Zooms computed from exp(p(7)*B2(:,:,1)+p(8)*B2(:,:,2)+p(9)*B2(:,:,3))
    %-----------------------------------------------------------------------
    B2        = zeros(4,4,3);
    B2(1,1,1) = 1;
    B2(2,2,2) = 1;
    B2(3,3,3) = 1;

    p      = zeros(9,1); % Parameters
    for it=1:10000,
        [R,dR] = spm_dexpm(p(1:6),B);  % Rotations + Translations
        [Z,dZ] = spm_dexpm(p(7:9),B2); % Zooms

        M  = R*Z; % Voxel-to-world estimate
        dM = zeros(4,4,6);
        for i=1:6, dM(:,:,i)   = dR(:,:,i)*Z; end
        for i=1:3, dM(:,:,i+6) = R*dZ(:,:,i); end
        dM = reshape(dM,[16,9]);

        d   = M(:)-M_avg(:); % Difference
        gr  = dM'*d;         % Gradient
        Hes = dM'*dM;        % Hessian
        p   = p - Hes\gr;    % Gauss-Newton update
        if sum(gr.^2)<1e-8, break; end
    end
    M_avg = M;
end

% Ensure that the FoV covers all images, with a few voxels to spare
%-----------------------------------------------------------------------
mn    =  Inf*ones(3,1);
mx    = -Inf*ones(3,1);
for i=1:size(Mat0,3),
    dm      = [dims(i,:) 1 1];
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1];
    M  = M_avg\Mat0(:,:,i);
    vx = M(1:3,:)*corners;
    mx = max(mx,max(vx,[],2));
    mn = min(mn,min(vx,[],2));
end
mx    = ceil(mx-(mx-mn)*0.05);
mn    = floor(mn+([mx(1:2)-mn(1:2);0])*0.05);
d     = (mx-mn+1)';
M_avg = M_avg * [eye(3) mn-1; 0 0 0 1];
M_avg(4,:)=[0 0 0 1];
return;
%==========================================================================

%==========================================================================
function B = se3_basis
% Basis functions for the lie algebra of the special Eucliden group
% (SE(3)).
B        = zeros(4,4,6);
B(1,4,1) = 1;
B(2,4,2) = 1;
B(3,4,3) = 1;
B([1,2],[1,2],4) = [0 1;-1 0];
B([3,1],[3,1],5) = [0 1;-1 0];
B([2,3],[2,3],6) = [0 1;-1 0];
return
%==========================================================================

%==========================================================================
function [y,E] = warp(x,M,phi,ord,varargin)
if ~isempty(varargin)
    d = varargin{1};
    if numel(varargin)>1
        st = varargin{2};
    else
        st = [1 1 1];
    end
    phi = Identity(d,st);
end
E = affine_transf(M,phi);

y = zeros(size(E(:,:,:,1)),'double');
for cl=1:size(x,4)
    x(:,:,:,cl) = spm_diffeo('bsplinc',x(:,:,:,cl),ord);
    y(:,:,:,cl) = spm_diffeo('bsplins',x(:,:,:,cl),E,ord);
end
y(~isfinite(y)) = 0;
%==========================================================================

%==========================================================================
function [Co,d1] = Identity(d,st,varargin)
d  = single(d);
st = single(st);
[e1,e2,e3] = meshgrid(st(2):st(2):d(2),st(1):st(1):d(1),st(3):st(3):d(3));
d1         = size(e1);
if ~isempty(varargin)
    Co      = cat(2,e2(:),e1(:),e3(:));
    Co(:,4) = 1;
else
    Co = cat(4,e2,e1,e3);
end
%==========================================================================