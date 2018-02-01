function [buf,nm,vr0] = init_buf(N,obj,V,x0,y0,z0,o,M,tpm,tot_S)
Kb       = numel(tpm.V);
d        = [size(x0) length(z0)];
modality = obj.modality;

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
else,
    error('Too many dimensions.');
end

nm = 0; % Number of voxels

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

cl   = cell(length(z0),1);
buf  = struct('msk',cl,'nm',cl,'f',cl,'dat',cl,'bf',cl,'code',cl);
for z=1:length(z0)
    if tot_S==1
        % Load only those voxels that are more than 5mm up
        % from the bottom of the tissue probability map.  This
        % assumes that the affine transformation is pretty close.
        z1  = M(3,1)*x0 + M(3,2)*y0 + (M(3,3)*z0(z) + M(3,4));
        e   = sqrt(sum(tpm.M(1:3,1:3).^2));
        e   = 5./e; % mm from edge of TPM
        buf(z).msk = z1>e(3);
    else
        buf(z).msk = ones(d(1:2),'logical');          
    end
    
    if isfield(obj,'msk') && ~isempty(obj.msk)
        % Exclude any voxels to be masked out
        msk        = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        buf(z).msk = buf(z).msk & msk; 
    end
    
    % Initially load all the data, but prepare to exclude
    % locations where any of the images is not finite, or
    % is zero.  We want this to work for skull-stripped
    % images too. The -3924 and -1500 options have been
    % added for CT data.
    fz = cell(1,N);
    for n=1:N
        fz{n}   = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        if strcmp(modality,'MRI'),    
            msk = isfinite(fz{n}) & (fz{n}~=0);
        elseif strcmp(modality,'CT'), 
            msk = isfinite(fz{n}) & (fz{n}~=min(fz{n}(:))) & (fz{n}~=0) & (fz{n}~=-3024) & (fz{n}~=-1500);            
        end
        buf(z).msk = buf(z).msk & msk;
    end
    buf(z).nm = nnz(buf(z).msk);
    nm        = nm + buf(z).nm; 
  
    code            = zeros([nnz(buf(z).msk) 1],typ);
    for n=1:N, code = bitor(code,bitshift(feval(cast,true([nnz(buf(z).msk) 1])),(n - 1))); end
    buf(z).code     = code;    

    % Eliminate unwanted voxels
    for n=1:N
        if scrand(n)
            % Data is an integer type, so to prevent aliasing in the histogram, small
            % random values are added.  It's not elegant, but the alternative would be
            % too slow for practical use.
            buf(z).f{n} = single(fz{n}(buf(z).msk) + rand(buf(z).nm,1)*scrand(n)-scrand(n)/2);
        else
            buf(z).f{n} = single(fz{n}(buf(z).msk{n}));
        end
        
        if strcmp(modality,'CT'), buf(z).f{n} = buf(z).f{n} + 1000; end
        
        mom0(n) = mom0(n) + buf(z).nm;
        mom1(n) = mom1(n) + sum(buf(z).f{n});
        mom2(n) = mom2(n) + sum(buf(z).f{n}.^2);
    end

    % Create a buffer for tissue probability info
    buf(z).dat = zeros([buf(z).nm,Kb],'single');
end

% Construct a ``Wishart-style prior'' (vr0)
vr0 = diag(mom2./mom0 - (mom1./mom0).^2)/Kb^2;
%=======================================================================