function [cls,M1] = write_res(obj)
lkp = obj.segment.lkp;
Kb  = max(lkp.part);
N   = numel(obj.image);

if ~isfield(obj.write_res,'write_tc'), obj.write_res.write_tc = true(Kb,4); end % native, import, warped, warped-mod
if ~isfield(obj.write_res,'write_bf'), obj.write_res.write_bf = false(N,2); end % field, corrected
if ~isfield(obj.write_res,'write_df'), obj.write_res.write_df = false(1,2); end % inverse, forward
if ~isfield(obj.write_res,'mrf'),      obj.write_res.mrf = 1;         end % MRF parameter
if ~isfield(obj.write_res,'cleanup_gwc'),  obj.write_res.cleanup_gwc = 1;     end % Run the ad hoc cleanup
if ~isfield(obj.write_res,'cleanup_lesion'),  obj.write_res.cleanup_lesion = 1;     end % Run the ad hoc cleanup
if ~isfield(obj.write_res,'bb'),       obj.write_res.bb = NaN(2,3);   end % Default to TPM bounding box
if ~isfield(obj.write_res,'vx'),       obj.write_res.vx = NaN;        end % Default to TPM voxel size
if ~isfield(obj.write_res,'G'),        obj.write_res.G = ones([Kb,1],'single'); end 
if ~isfield(obj.write_res,'nmrf_its'), obj.write_res.nmrf_its = 10; end 
if ~isfield(obj,'dir_write'),          obj.dir_write = [];       end % Output directory

write_tc = obj.write_res.write_tc;
write_bf = obj.write_res.write_bf;
write_df = obj.write_res.write_df;
mrf      = obj.write_res.mrf;
cleanup_gwc  = obj.write_res.cleanup_gwc;
cleanup_lesion  = obj.write_res.cleanup_lesion;
bb       = obj.write_res.bb;
vx       = obj.write_res.vx;
G        = obj.write_res.G*mrf;
nmrf_its = obj.write_res.nmrf_its;
odir     = obj.dir_write;
modality = obj.modality;

% Get parameters
%--------------------------------------------------------------------------
wp_l = obj.segment.wp_l;

% Load template
%-----------------------------------------------------------------------
tpm = spm_load_logpriors(obj.pth_template,obj.segment.wp);

% Read essentials from tpm (it will be cleared later)
d1 = size(tpm.dat{1});
if numel(d1)==2
    d1(3) = 1;
end
M1 = tpm.M;

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


[pth,nam] = fileparts(obj.image(1).fname);
if ~isempty(odir) && ischar(odir), pth = odir; end
ind  = obj.image(1).n;
d    = obj.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3 = 1:d(3);

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N
    d3         = [size(obj.segment.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = obj.segment.Tbias{n};

    % Need to fix writing of bias fields or bias corrected images, when the data used are 4D.
    [pth1,nam1] = fileparts(obj.image(n).fname);
    if ~isempty(odir) && ischar(odir), pth1 = odir; end
    chan(n).ind = obj.image(n).n;

    if write_bf(n,2)
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth1,['m', nam1, '.nii']),...
                                     obj.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nc.mat  = obj.image(n).mat;
        chan(n).Nc.mat0 = obj.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end

    if write_bf(n,1),
        chan(n).Nf      = nifti;
        chan(n).Nf.dat  = file_array(fullfile(pth1,['BiasField_', nam1, '.nii']),...
                                     obj.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nf.mat  = obj.image(n).mat;
        chan(n).Nf.mat0 = obj.image(n).mat;
        chan(n).Nf.descrip = 'Estimated Bias Field';
        create(chan(n).Nf);
    end
end

do_cls   = any(write_tc(:)) || nargout>=1;
tiss(Kb) = struct('Nt',[]);
for k1=1:Kb
    if write_tc(k1,4) || any(write_tc(:,3)) || write_tc(k1,2) || nargout>=1,
        do_cls  = true;
    end
    if write_tc(k1,1),
        tiss(k1).Nt      = nifti;
        tiss(k1).Nt.dat  = file_array(fullfile(pth,['c', num2str(k1), nam, '.nii']),...
                                      obj.image(n).dim(1:3),...
                                      [spm_type('uint8') spm_platform('bigend')],...
                                      0,1/255,0);
        tiss(k1).Nt.mat  = obj.image(n).mat;
        tiss(k1).Nt.mat0 = obj.image(n).mat;
        tiss(k1).Nt.descrip = ['Tissue class ' num2str(k1)];
        create(tiss(k1).Nt);
        do_cls = true;
    end
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);

Nii   = nifti(obj.pth_vel);
Twarp = single(Nii.dat(:,:,:,:)); 

Coef{1} = spm_bsplinc(Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(Twarp(:,:,:,3),prm);

do_defs = any(write_df);
do_defs = do_cls | do_defs;
if do_defs
    if write_df(1)
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam, '.nii']),...
                               [obj.image(1).dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
        Ndef.mat  = obj.image(1).mat;
        Ndef.mat0 = obj.image(1).mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
    if write_df(2) || any(any(write_tc(:,[2,3,4]))) || nargout>=1
        y = zeros([obj.image(1).dim(1:3),3],'single');
    end
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\obj.Affine*obj.image(1).mat;

%--------------------------------------------------------------------------
if obj.segment.mrf.do_mrf  
    resp = init_resp(obj,lkp,d);
    
    % Compute neighborhood part of responsibilities
    % (from previous responsibilities)
    lnPzN = compute_lnPzN(obj.segment.mrf,resp);
    
    clear resp
end

%--------------------------------------------------------------------------
if do_cls
    Q = zeros([d(1:3),Kb],'single');
end

for z=1:length(x3)

    % Bias corrected image
    f   = cell(1,N);
    bf  = cell(1,N);
    msk = cell(1,N);
    for n=1:N
        f{n}   = spm_sample_vol(obj.image(n),x1,x2,o*x3(z),0);        
        msk{n} = spm_misc('msk_modality',f{n},obj.modality);        
        
        f{n}(~msk{n}) = 0;
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
        bfn = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        
        if ~isempty(chan(n).Nc),
            % Write a plane of bias corrected data
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bfn.*f{n};
        end
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bfn;
        end
        
%         f{n}  = f{n}(msk{n});        
        bf{n} = bfn(:);             
    end    
    clear bfn
    
%     % Bias corrected image
%     f   = cell(1,N);
%     bf  = cell(1,N);
%     msk = cell(1,N);
%     for n=1:N
%         f{n}   = spm_sample_vol(obj.image(n),x1,x2,o*x3(z),0);
%         msk{n} = spm_misc('msk_modality',f{n},obj.modality);              
%         bf{n}  = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));  
%         
%         if ~isempty(chan(n).Nc),
%             % Write a plane of bias corrected data
%             chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf{n}.*f{n};
%         end
%         if ~isempty(chan(n).Nf),
%             % Write a plane of bias field
%             chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf{n};
%         end
%     end    
%         
%     if ~obj.segment.do_missing_data
%         tmp = true;
%         for n=1:N
%             tmp = tmp & msk{n};
%         end
% 
%         for n=1:N
%             msk{n} = tmp;
%         end
%     end
    
    if do_defs
        % Compute the deformation (mapping voxels in image to voxels in TPM)
        [t1,t2,t3] = make_deformation(Coef,z,obj.segment.MT,prm,x1,x2,x3,M);

        if exist('Ndef','var')
            % Write out the deformation to file, adjusting it so mapping is
            % to voxels (voxels in image to mm in TPM)
            Ndef.dat(:,:,z,1,1) = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,2) = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,3) = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
        end

        if exist('y','var')
            % If needed later, save in variable y
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        if do_cls
%             code            = zeros([numel(msk{1}) 1],typ);
%             for n=1:N, code = bitor(code,bitshift(feval(cast,msk{n}(:)),(n - 1))); end                   
% 
%             % Parametric representation of intensity distributions                       
%             q  = zeros([d(1:2) Kb]);                          
%             q1 = log_likelihoods(f,bf,obj.segment.mg,obj.segment.gmm,msk,code);
%             q1 = reshape(q1,[d(1:2),numel(obj.segment.mg)]);
%             q1 = exp(q1) + eps;
%             
%             b = spm_sample_logpriors(tpm,t1,t2,t3);    
%             s = zeros(size(b{1}));
%             for k1 = 1:Kb
%                 b{k1} = obj.segment.wp(k1)*b{k1};
%                 s     = s + b{k1};
%             end
%             
%             for k1=1:Kb
%                 tmp                 = sum(q1(:,:,lkp.part==k1),3);
%                 tmp(~isfinite(tmp)) = 1e-3;
%                 q(:,:,k1)           = tmp.*(b{k1}./s);
%             end
%                 
%             Q(:,:,z,:) = reshape(q,[d(1:2),1,Kb]);

            code            = zeros([numel(msk{1}) 1],typ);
            for n=1:N, code = bitor(code,bitshift(feval(cast,msk{n}(:)),(n - 1))); end  

            msk1 = code>0;
            b    = spm_sample_logpriors(tpm,t1(:),t2(:),t3(:));
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

            B = zeros([numel(msk1) Kb]);
            for k1=1:Kb, 
                B(:,k1) = b{k1}(:); 
            end
            clear b msk1

            % Get neighborhood term
            if obj.segment.mrf.do_mrf  
                lnPzN1 = double(reshape(lnPzN(:,:,z,:),[prod(d(1:2)) Kb]));
            else
                lnPzN1 = zeros([1 Kb]);
            end

            q1 = latent(f,bf,obj.segment.mg,obj.segment.gmm,B,lnPzN1,lkp,obj.segment.wp,code,labels,wp_l);
            clear B f bf lnPzN1

            q = zeros([prod(d(1:2)) Kb]);  
            for k1=1:Kb
                tmp                 = sum(q1(:,lkp.part==k1),2);
                tmp(~isfinite(tmp)) = 1e-3;
                q(:,k1)             = tmp;
            end 
            clear q1 tmp

            Q(:,:,z,:) = single(reshape(q,[d(1:2),1,Kb]));
            clear q
        end
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');
clear lnPzN

cls   = cell(1,Kb);
if do_cls
    P = zeros([d(1:3),Kb],'uint8');
    if mrf==0
        % Normalise to sum to 1
        sQ = (sum(Q,4)+eps)/255;
        for k1=1:size(Q,4)
            P(:,:,:,k1) = uint8(round(Q(:,:,:,k1)./sQ));
        end
        clear sQ
    else
        % Use a MRF cleanup procedure        
        spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
        vx2 = 1./single(sum(obj.image(1).mat(1:3,1:3).^2));
        for iter=1:nmrf_its
            spm_mrf(P,Q,G,vx2);
            spm_progress_bar('set',iter);
        end
    end
    clear Q

    if cleanup_gwc
        % Use an ad hoc brain cleanup procedure
        if size(P,4)>3
            P = clean_gwc(P,cleanup_gwc);
        else
            warning('Cleanup not done.');
        end
    end

    if cleanup_lesion
        % Use an ad hoc lesion cleanup procedure
        if size(P,4)>3
            P = clean_lesion(P,cleanup_lesion);
        end
    end
    
    % Write tissues if necessary
    for k1=1:Kb
        if ~isempty(tiss(k1).Nt)
            for z=1:length(x3)
                tmp = double(P(:,:,z,k1))/255;
                tiss(k1).Nt.dat(:,:,z,ind(1),ind(2)) = tmp;
            end
        end
    end
    spm_progress_bar('clear');

    % Put tissue classes into a cell array...
    for k1=1:Kb
        if write_tc(k1,4) || any(write_tc(:,3)) || write_tc(k1,2) || nargout>=1
            cls{k1} = P(:,:,:,k1);
        end
    end
    clear P % ...and remove the original 4D array
end

clear tpm
M0  = obj.image(1).mat;

if any(write_tc(:,2))
    % "Imported" tissue class images

    % Generate mm coordinates of where deformations map from
    x      = affind(rgrid(d),M0);

    % Generate mm coordinates of where deformation maps to
    y1     = affind(y,M1);

    % Procrustes analysis to compute the closest rigid-body
    % transformation to the deformation, weighted by the
    % interesting tissue classes.
    ind        = find(write_tc(:,2)); % Saved tissue classes
    [dummy,R]  = spm_get_closest_affine(x,y1,single(cls{ind(1)})/255);
    clear x y1

    mat0   = R\mat; % Voxel-to-world of original image space

    fwhm   = max(vx./sqrt(sum(obj.image(1).mat(1:3,1:3).^2))-1,0.01);
    for k1=1:size(write_tc,1)
        if write_tc(k1,2)

            % Low pass filtering to reduce aliasing effects in downsampled images,
            % then reslice and write to disk
            tmp1    = decimate(single(cls{k1}),fwhm);
            Ni      = nifti;
            Ni.dat  = file_array(fullfile(pth,['rc', num2str(k1), nam, '.nii']),...
                                 odim,...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
            Ni.mat         = mat;
            Ni.mat_intent  = 'Aligned';
            Ni.mat0        = mat0;
            Ni.mat0_intent = 'Aligned';
            Ni.descrip     = ['Imported Tissue ' num2str(k1)];
            create(Ni);

            for i=1:odim(3)
                tmp           = spm_slice_vol(tmp1,M0\mat0*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
                Ni.dat(:,:,i) = tmp;
            end
            clear tmp1
        end
    end
end

if any(write_tc(:,3)) || any(write_tc(:,4)) || nargout>=1 || write_df(2)
    % Adjust stuff so that warped data (and deformations) have the
    % desired bounding box and voxel sizes, instead of being the same
    % as those of the tissue probability maps.
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
end

if any(write_tc(:,3)) || any(write_tc(:,4)) || nargout>=1

    if any(write_tc(:,3))
        C = zeros([d1,Kb],'single');
    end

    spm_progress_bar('init',Kb,'Warped Tissue Classes','Classes completed');
    for k1 = 1:Kb
        if ~isempty(cls{k1})
            c = single(cls{k1})/255;
            if any(write_tc(:,3))
                [c,w]       = spm_diffeo('push',c,y,d1(1:3));
                vx          = sqrt(sum(M1(1:3,1:3).^2));
                spm_field('boundary',1);
                C(:,:,:,k1) = spm_field(w,c,[vx  1e-6 1e-4 0  3 2]);
                clear w
            else
                c      = spm_diffeo('push',c,y,d1(1:3));
            end
            if nargout>=1
                cls{k1} = c;
            end
            if write_tc(k1,4)
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['mwc', num2str(k1), nam, '.nii']),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = c*abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
            end
            spm_progress_bar('set',k1);
        end
    end
    spm_progress_bar('Clear');

    if any(write_tc(:,3))
        spm_progress_bar('init',Kb,'Writing Warped Tis Cls','Classes completed');
        C = max(C,eps);
        s = sum(C,4);
        for k1=1:Kb
            if write_tc(k1,3)
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['wc', num2str(k1), nam, '.nii']),...
                                    d1,'uint8',0,1/255,0);
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = C(:,:,:,k1)./s;
            end
            spm_progress_bar('set',k1);
        end
        spm_progress_bar('Clear');
        clear C s
    end    
end

if write_df(2)
    y         = spm_diffeo('invdef',y,d1,eye(4),M0);
    y         = spm_extrapolate_def(y,M1);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam, '.nii']),...
                           [d1,1,3],'float32',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

return;
%==========================================================================

%==========================================================================
% function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
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
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
return;
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
return;
%==========================================================================

%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
return;
%==========================================================================

%=======================================================================
function R = latent(f,bf,mg,gmm,lnPa,lnPzN,lkp,wp,code,labels,wp_l,cr)
if nargin<12, cr = []; end

Kb   = max(lkp.part); 
tiny = 1e-4; 

lnPa = log_spatial_priors(lnPa,wp);
lnPx = log_likelihoods(f,bf,mg,gmm,code,cr);

if ~isempty(labels) 
    lnPa1 = zeros([size(lnPa,1) Kb]); 
end

msk = code>0;
R   = zeros(size(lnPa));
for k1=1:Kb
    for k=find(lkp.part==k1)
        if ~isempty(labels)
            if lkp.lab(k1)~=0               
                msk_l = labels==lkp.lab(k1);
 
                beta1 = log(1 - tiny); % log-probability of true labels
                beta2 = log(tiny);    % log-probability of false labels
 
                lnPa1(msk_l,k1)  = (1 - wp_l)*lnPa(msk_l,k1) + wp_l*beta1;
                lnPa1(~msk_l,k1) = (1 - wp_l)*lnPa(~msk_l,k1) + wp_l*beta2;                
            else
                lnPa1(:,k1) = lnPa(:,k1);
            end
            
            if numel(lnPzN)==Kb
                R(:,k)   = lnPa1(:,k1) + lnPzN(:,k1);
                R(msk,k) = R(msk,k) + lnPx(msk,k);
            else
                R(:,k) = lnPa1(:,k1) + lnPzN(:,k1);
                R(msk,k) = R(msk,k) + lnPx(msk,k);
            end            
        else
            if numel(lnPzN)==Kb
                R(:,k) = lnPa(:,k1) + lnPzN(:,k1);
                R(msk,k) = R(msk,k) + lnPx(msk,k);
            else
                R(:,k) = lnPa(:,k1) + lnPzN(:,k1);
                R(msk,k) = R(msk,k) + lnPx(msk,k);
            end
        end
    end
end

if ~isempty(labels) 
    lnPa = lnPa1;
    clear B1
end

logSumQ = spm_matcomp('logsumexp',R,2);
logQ    = bsxfun(@minus,R,logSumQ);
R       = exp(logQ);
%=======================================================================

%=======================================================================
function L = log_likelihoods(f,bf,mg,gmm,code,cr)
if nargin<6, cr = []; end

K = numel(mg);
N = numel(f);
I = numel(f{1});

% Compute Bx
%--------------------------------------------------------------------------
if isempty(cr)
    cr                       = zeros(I,N);
    for n=1:N, cr(:,n)  = double(f{n}(:)).*double(bf{n}(:)); end
elseif iscell(cr)
    cr1                      = zeros(I,N);
    for n=1:N, cr1(:,n) = cr{n}; end    
    cr                       = cr1; clear cr1
end

% Compute log|B|
%--------------------------------------------------------------------------
nbf                      = zeros([I N]);
for n=1:N, nbf(:,n) = double(bf{n}); end
clear bf

% Compute likelihoods
%--------------------------------------------------------------------------
L = zeros(I,K);
for n=2:2^N
    msk0 = dec2bin(n - 1,N)=='1';
    ind  = find(code==msk0*(2.^(0:(N - 1))'));
    if ~isempty(ind)
        for k=1:K            
            d        = bsxfun(@minus,cr(ind,msk0)',gmm.po.m(msk0,k));
            Q        = chol(gmm.po.W(msk0,msk0,k))*d;
            E        = N/gmm.po.b(k) + gmm.po.n(k)*dot(Q,Q,1);
            L(ind,k) = 0.5*(spm_prob('Wishart','Elogdet',gmm.po.W(msk0,msk0,k),gmm.po.n(k)) - E') + log(mg(k)) - N/2*log(2*pi) + log(prod(nbf(ind,msk0),2));            
        end
    end
end
%=======================================================================