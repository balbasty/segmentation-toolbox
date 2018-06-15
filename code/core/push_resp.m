function obj = push_resp(obj)
lkp      = obj.segment.lkp;
Kb       = max(lkp.part);
N        = numel(obj.image);
modality = obj.modality;
mrf      = obj.segment.mrf;

% Get parameters
%--------------------------------------------------------------------------
wp_l = obj.segment.wp_l;

% Load template
%--------------------------------------------------------------------------
tpm = spm_load_logpriors(obj.pth_template,obj.segment.wp);
M1  = tpm.M;

% For missing data
%--------------------------------------------------------------------------
if N<=8,
    cast = @uint8;
    typ  = 'uint8';
elseif N<=16,
    cast = @uint16;
    typ  = 'uint16';
elseif N<=32,
    cast = @uint32;
    typ  = 'uint32';
elseif N<=64,
    cast = @uint64;
    typ  = 'uint64';
else,
    error('Too many dimensions.');
end

% Define orientation and field of view of any "normalised" space
% data that may be generated (wc*.nii, mwc*.nii, rc*.nii & y_*.nii).
% Use the actual dimensions and orientations of the tissue priors.
%--------------------------------------------------------------------------
odim = tpm.V(1).dim;
mat  = tpm.V(1).mat;

%--------------------------------------------------------------------------
d         = obj.image(1).dim(1:3);
[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3        = 1:d(3);

%--------------------------------------------------------------------------
chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N
    d3         = [size(obj.segment.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = obj.segment.Tbias{n};

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
if mrf.do_mrf        
    resp = init_resp(obj,lkp,d);
    
    % Compute neighborhood part of responsibilities
    % (from previous responsibilities)
    lnPzN = compute_lnPzN(mrf,resp);
    
    clear resp
end

%--------------------------------------------------------------------------
y = zeros([obj.image(1).dim(1:3),3],'single');
Q = NaN([d(1:3),Kb],'single');
for z=1:length(x3)
    % Bias corrected image
    f   = cell(1,N);
    bf  = cell(1,N);
    nm  = 0;    
    msk = cell(1,N);
    for n=1:N
        f{n}   = spm_sample_vol(obj.image(n),x1,x2,o*x3(z),0);        
        msk{n} = spm_misc('msk_modality',f{n},obj.modality);        
    end
    
    if ~obj.segment.do_missing_data
        tmp = true;
        for n=1:N
            tmp = tmp & msk{n};
        end

        for n=1:N
            msk{n} = tmp;
        end
    end
    
    for n=1:N
        f{n}  = f{n}(msk{n});
        bfn   = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bf{n} = bfn(msk{n});      
        
        nm = nm + nnz(msk{n});
    end    
    clear bfn
       
    % Compute the deformation (mapping voxels in image to voxels in TPM)
    [t1,t2,t3] = make_inv_deformation(Coef,z,obj.segment.MT,prm,x1,x2,x3,M);

    y(:,:,z,1) = t1;
    y(:,:,z,2) = t2;
    y(:,:,z,3) = t3;    
    
    if nm==0
        continue;
    end
        
    code            = zeros([numel(msk{1}) 1],typ);
    for n=1:N, code = bitor(code,bitshift(feval(cast,msk{n}(:)),(n - 1))); end  
        
    msk1 = code>0;
    b    = spm_sample_logpriors(tpm,t1(msk1),t2(msk1),t3(msk1));
    clear t1 t2 t3
    
    if isempty(lkp.lab) || isempty(obj.labels)
        labels = [];
    else
        tmp = uint8(obj.labels.private.dat(:,:,z));        
        tmp = tmp(msk1);        
        tmp = tmp(:); 
        
        msk2       = ismember(tmp,lkp.lab);
        tmp(~msk2) = 0;
        
        labels = zeros([numel(tmp) 1],'uint8');         
        for k=1:Kb      
            msk2         = tmp==lkp.lab(k) & lkp.lab(k)~=0; 
            labels(msk2) = obj.segment.lkp.lab(k); 
        end 
        clear tmp 
    end
    
    B = zeros([nnz(msk1) Kb]);
    for k1=1:Kb, 
        B(:,k1) = b{k1}(:); 
    end
    clear b msk1
    
    % Get neighborhood term
    if mrf.do_mrf   
        lnPzN1 = double(reshape(lnPzN(:,:,z,:),[prod(d(1:2)) Kb]));
    else
        lnPzN1 = zeros([1 Kb]);
    end
    
    q1 = latent(f,bf,obj.segment.mg,obj.segment.gmm,B,lnPzN1,lkp,obj.segment.wp,msk,code,labels,wp_l);
    clear B f bf lnPzN1
    
    q = NaN([prod(d(1:2)) Kb]);  
    for k1=1:Kb
        q(:,k1) = sum(q1(:,lkp.part==k1),2);
    end 
    clear q1 
    
    Q(:,:,z,:) = single(reshape(q,[d(1:2),1,Kb]));
    clear q
end
clear tpm x1 x2 x3 o Coeff chan msk lnPzN

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
for k1=1:Kb
    [t1,c] = spm_diffeo('push',Q(:,:,:,k1),y,d1(1:3));
     
    if k1==1
        % Calculate bounding box from count image
        bb_push          = compute_bb(c,obj.pth_resp{k1},M1);        
        rngx             = bb_push(1,1):bb_push(1,2);
        rngy             = bb_push(2,1):bb_push(2,2);
        rngz             = bb_push(3,1):bb_push(3,2);
        dm_bb            = [numel(rngx) numel(rngy) numel(rngz)];
        obj.push_resp.bb = bb_push;
        
        % Allocate sufficient statistics
        t = zeros([dm_bb Kb],'single');
    end
    
    fname = obj.pth_resp{k1};
    if (exist(fname,'file')==2), delete(fname); end
    
    Nii         = nifti;
    Nii.dat     = file_array(fname,dm_bb,obj.dt,0,1,0);                    
    Nii.mat     = M1;    
    Nii.mat0    = M1;    
    Nii.descrip = ['Pushed responsibilities ' num2str(k1)];        
    create(Nii);
    
    t(:,:,:,k1)    = t1(rngx,rngy,rngz);    
    Nii.dat(:,:,:) = t(:,:,:,k1);    
end
clear Q y Nii t1 c

% Prepare computing template gradient and Hessian
%--------------------------------------------------------------------------
Nii  = nifti(obj.pth_template);
d    = [dm_bb,Kb,1,1,1];
R    = null(ones(1,d(4)));    
lwp  = reshape(log(obj.segment.wp),1,1,d(4));
lwp1 = reshape(log(obj.segment.wp),1,1,1,d(4));

t(:,:,:,lkp.keep) = max(t(:,:,:,lkp.keep),eps('single')*1000);
s                 = sum(t,4);

% Compute gradients and Hessian
%--------------------------------------------------------------------------
W   = zeros([d(1:3) round(((d(4) - 1)*d(4))/2)],'single'); % 2nd derivatives
gr  = zeros([d(1:3),d(4) - 1],'single');    

ll  = 0;
for z=1:d(3), % Loop over planes

    % Log of template    
    sz = Nii.dat(rngx,rngy,rngz(z),:);
    
    % log-likelihood (using log-sum-exp)
    sm0 = bsxfun(@plus,sz,lwp1);            
    
    t1 = t(:,:,z,:);    
    if ~isempty(lkp.lab) 
        for k1=1:Kb 
            if lkp.lab(k1)~=0                
               t1(:,:,1,k1) = (1 - wp_l)*t1(:,:,1,k1);   
            end 
        end         
    end 
        
    sm1 = sum(t1.*sm0,4);
    
    t1 = sum(t1,4);
    
    sm2 = spm_matcomp('logsumexp',sm0,4).*t1;
    
    ll0 = sum(sum(sm1 - sm2));                        
    ll  = ll - ll0;
    
    sz = reshape(sz,d(1),d(2),d(4));
    a  = zeros([d(1:2) 1 d(4) - 1],'single');
    for j1=1:(d(4)-1),
        az = zeros(d(1:2));
        for j2=1:d(4),
            az = az + R(j2,j1)*sz(:,:,j2); % Note the rotation
        end
        a(:,:,:,j1) = az;
    end

    % Compute softmax for this plane
    mu = double(reshape(sftmax(a,R,lwp),[d(1:2),d(4)]));
    
%     if ~isempty(lkp.lab) 
%         for k1=1:Kb 
%             if lkp.lab(k1)~=0                
%                mu(:,:,k1) = mu(:,:,k1).^(1 - wp_l);
%             end 
%         end         
%     end    

    % Compute first derivatives (d(4)-1) x 1 
    grz = bsxfun(@times,mu,s(:,:,z)) - double(reshape(t(:,:,z,:),[d(1:2),d(4)]));
    
    if ~isempty(lkp.lab) 
        for k1=1:Kb 
            if lkp.lab(k1)~=0                
               grz(:,:,k1) = grz(:,:,k1)*(1 - wp_l);
            end 
        end         
    end    

    for j1=1:(d(4)-1),
        gr1 = zeros([d(1:2) 1 d(4) - 1],'single');
        for j2=1:d(4),
            gr1(:,:,1,j1) = gr1(:,:,1,j1) + R(j2,j1)*grz(:,:,j2); % Note the rotation
        end
        gr(:,:,z,j1) = gr(:,:,z,j1) + gr1(:,:,1,j1);%.*s(:,:,z);
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

    if ~isempty(lkp.lab) 
        for k1=1:Kb 
            for k2=1:Kb
                if lkp.lab(k1)~=0 || lkp.lab(k2)~=0
                   wz(:,:,k1,k2) = wz(:,:,k1,k2)*(1 - wp_l);
                end 
            end
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
        W(:,:,z,j1) = W(:,:,z,j1) + wz(:,:,j1,j1) + 1e-4*d(4)^2;
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
clear t s mu wz1 a Nii

obj.ll_template = ll;  

% Save gradients
%--------------------------------------------------------------------------
d = size(gr);
for k1=1:d(4)
    fname = obj.pth_gr{k1};
    if (exist(fname,'file')==2), delete(fname); end
    
    Nii      = nifti;
    Nii.dat  = file_array(fname,d(1:3),obj.dt,0,1,0);
    Nii.mat  = M1;
    Nii.mat0 = M1;
    Nii.descrip = ['Template gradients ' num2str(k1)];
    create(Nii);
    
    Nii.dat(:,:,:) = gr(:,:,:,k1);   
end
clear Nii Nii1 gr

% Save Hessian
%--------------------------------------------------------------------------
d = size(W);
for k1=1:d(4)
    fname = obj.pth_H{k1};
    if (exist(fname,'file')==2), delete(fname); end
    
    Nii      = nifti;
    Nii.dat  = file_array(fname,d(1:3),obj.dt,0,1,0);
    Nii.mat  = M1;
    Nii.mat0 = M1;
    Nii.descrip = ['Template Hessian ' num2str(k1)];
    create(Nii);
    
    Nii.dat(:,:,:) = W(:,:,:,k1);        
end
clear Nii Nii1 W

if ~isempty(obj.dir_der_local)
    % SSH copy from Holly to local machine
    cmd = ['scp -r ' fullfile(obj.dir_der,'*') ' ' obj.address_local ':' obj.dir_der_local];
    system(cmd);
end
%==========================================================================

%==========================================================================
function [x1,y1,z1] = make_inv_deformation(sol,z,MT,prm,x0,y0,z0,M)
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
if numel(z0)==1
   z1 = ones(size(z1));
end
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
function bb = compute_bb(img,pth,M)
[pth,nam,ext] = fileparts(pth);
nfname        = fullfile(pth,['temp-' nam ext]);
   
Nii         = nifti;
Nii.dat     = file_array(nfname,size(img),'float32',0,1,0);                    
Nii.mat     = M;    
Nii.mat0    = M;    
Nii.descrip = 'temp';        
create(Nii);

Nii.dat(:,:,:) = img;

V        = spm_vol(nfname);
[bb,vx1] = spm_get_bbox(V,'nz');

% vx   = abs(prod(vx1))^(1/3);
% odim = abs(round((bb_push(2,1:3)-bb_push(1,1:3))/vx)) + 1

% bb_push(1,:) = vx*round(bb_push(1,:)/vx);
% bb_push(2,:) = vx*round(bb_push(2,:)/vx);
    
mn = [bb(1,:) 1]';
mx = [bb(2,:) 1]';

mn = M\mn;
mx = M\mx;

bb = [mn(1:3), mx(1:3)];
bb = sort(bb,2);

delete(nfname); 
%==========================================================================