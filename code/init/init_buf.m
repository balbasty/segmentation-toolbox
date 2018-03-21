function [buf,nm,vr0,mn,mx,scl_int] = init_buf(N,obj,V,x0,y0,z0,o,M,tpm,tot_S)
Kb              = max(obj.lkp.part);
d               = [size(x0) length(z0)];
modality        = obj.modality;
do_missing_data = obj.do_missing_data;
uniform         = obj.uniform;

if isfield(obj,'msk') && ~isempty(obj.msk)
    VM = spm_vol(obj.msk);
    if sum(sum((VM.mat-V(1).mat).^2)) > 1e-6 || any(VM.dim(1:3) ~= V(1).dim(1:3))
        error('Mask must have the same dimensions and orientation as the image.');
    end
end

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
else
    error('Too many dimensions.');
end

nm = 0; % Number of voxels
mn = Inf*ones(N,1);
mx = -Inf*ones(N,1);

scrand = zeros(N,1);
for n=1:N
    if spm_type(V(n).dt(1),'intt')
        scrand(n) = V(n).pinfo(1);
    end
end

% Overall moments used later for regularising via a ``Wishart-style prior''
mom0 = zeros(1,N);
mom1 = zeros(1,N);
mom2 = zeros(1,N);  

sint    = zeros(1,N);
nms     = zeros(1,N);
scl_int = zeros(1,N);

cl   = cell(length(z0),1);
buf  = struct('msk',cl,'nm',cl,'Nm',cl,'f',cl,'dat',cl,'bf',cl,'code',cl,'img',cl,'labels',cl);
for z=1:length(z0)
    if tot_S==1
        % Load only those voxels that are more than 5mm up
        % from the bottom of the tissue probability map.  This
        % assumes that the affine transformation is pretty close.
        z1            = M(3,1)*x0 + M(3,2)*y0 + (M(3,3)*z0(z) + M(3,4));
        e             = sqrt(sum(tpm.M(1:3,1:3).^2));
        e             = 5./e; % mm from edge of TPM
        for n=1:N
            buf(z).msk{n} = z1>e(3);
        end
    else
        for n=1:N
            buf(z).msk{n} = ones(d(1:2),'logical');          
        end
    end
    
    if isfield(obj,'msk') && ~isempty(obj.msk)
        % Exclude any voxels to be masked out
        msk = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        for n=1:N
            buf(z).msk{n} = buf(z).msk{n} & msk; 
        end
    end
    
    % Load the data
    fz  = cell(1,N);
    msk = true;
    for n=1:N
        fz{n}         = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        buf(z).msk{n} = msk_modality(fz{n},modality,obj.trunc_ct);
        buf(z).nm(n)  = nnz(buf(z).msk{n});
                
        if uniform
            buf(z).img{n} = single(V(n).private.dat(:,:,z));
            msk           = msk & msk_modality(buf(z).img{n},modality,obj.trunc_ct);
        end
    end              
    
    if uniform
        for n=1:N            
            buf(z).img{n} = buf(z).img{n}(msk);
        end
    end
    
    if ~do_missing_data
        msk = true;
        for n=1:N
            msk = msk & buf(z).msk{n};
        end

        for n=1:N
            buf(z).msk{n} = msk;
            buf(z).nm(n)  = nnz(buf(z).msk{n});
        end
    end    
    
    code            = zeros([numel(buf(z).msk{1}) 1],typ);
    for n=1:N, code = bitor(code,bitshift(feval(cast,buf(z).msk{n}(:)),(n - 1))); end
    buf(z).code     = code;    

    buf(z).Nm = nnz(code);
    nm        = nm + buf(z).Nm; 
        
    % Prepare labels (if provided)
    %----------------------------------------------------------------------
    if isempty(obj.lkp.lab)
        buf(z).labels = [];
    else
        msk = code>0;
        tmp = uint8(spm_sample_vol(obj.labels,x0,y0,o*z0(z),0));        
        tmp = tmp(msk);
        clear msk
        
        nlabels = zeros([numel(tmp) Kb],'uint8');
        for k=1:Kb
            nlabels(:,k) = tmp==obj.lkp.lab(k) & obj.lkp.lab(k)~=0;             
        end
        clear tmp
        
        buf(z).labels = nlabels;
        clear nlabels
    end
    
    % Eliminate unwanted voxels
    %----------------------------------------------------------------------
    for n=1:N
        if scrand(n)
            % Data is an integer type, so to prevent aliasing in the histogram, small
            % random values are added.  It's not elegant, but the alternative would be
            % too slow for practical use.
            buf(z).f{n} = single(fz{n}(buf(z).msk{n}) + rand(buf(z).nm(n),1)*scrand(n)-scrand(n)/2);
            
            if uniform
                buf(z).img{n} = buf(z).img{n} + rand(size(buf(z).img{n}))*scrand(n)-scrand(n)/2;
            end
        else
            buf(z).f{n} = single(fz{n}(buf(z).msk{n}));
        end
        
        if max(buf(z).f{n}(:))>mx(n)
            mx(n) = max(buf(z).f{n}(:));
        end
        
        if min(buf(z).f{n}(:))<mn(n)
            mn(n) = min(buf(z).f{n}(:));
        end
        
        mom0(n) = mom0(n) + buf(z).nm(n);
        mom1(n) = mom1(n) + sum(buf(z).f{n});
        mom2(n) = mom2(n) + sum(buf(z).f{n}.^2);
        
        sint(n) = sint(n) + sum(buf(z).f{n});
        nms(n)  = nms(n)  + buf(z).nm(n);
    end

    % Create a buffer for tissue probability info
    buf(z).dat = zeros([buf(z).Nm,Kb],'single');
end

for n=1:N
    scl_int(n) = (1024 / (sint(n)/nms(n)));
end

% Construct a ``Wishart-style prior'' (vr0)
vr0 = diag(mom2./mom0 - (mom1./mom0).^2)/Kb^2;
%=======================================================================