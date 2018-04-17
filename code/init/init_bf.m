function varargout = init_bf(varargin)
buf = varargin{1};
obj = varargin{2};
V   = varargin{3};
x0  = varargin{4};
y0  = varargin{5};
z0  = varargin{6};
ff  = varargin{7};
scl = varargin{8};

N    = numel(buf(1).f);
cl   = cell(N,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});
for n=1:N
    biasreg = obj.segment.biasreg(n)*ff; 
    vx      = spm_misc('vxsize',V(n).mat);
    fwhm    = obj.segment.biasfwhm(n);
    d0      = V(n).dim;
    
    sd = vx(1)*d0(1)/fwhm; d3(1) = ceil(sd*2); krn_x = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
    sd = vx(2)*d0(2)/fwhm; d3(2) = ceil(sd*2); krn_y = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
    sd = vx(3)*d0(3)/fwhm; d3(3) = ceil(sd*2); krn_z = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
    
    % BENDING ENERGY regularisation for bias correction (when building TPMs)
    chan(n).C = spm_sparse('precision','field',[numel(krn_x) numel(krn_y) numel(krn_z)],vx,[0 0 biasreg 0 0],'n');
    
    % Initial parameterisation of bias field
    if isfield(obj.segment,'Tbias') && ~isempty(obj.segment.Tbias{n})                
        chan(n).T = obj.segment.Tbias{n};
        
        if obj.tot_S>1
            % When building TPMs, mean correct DC component of bias field
            chan(n).T(1,1,1) = chan(n).T(1,1,1) - obj.segment.bf.avg_dc(n); 
        end
    else
        chan(n).T = zeros(d3);        
    end        
    
    % Basis functions for bias correction
    chan(n).B3 = spm_dctmtx(d0(3),d3(3),z0);
    chan(n).B2 = spm_dctmtx(d0(2),d3(2),y0(1,:)');
    chan(n).B1 = spm_dctmtx(d0(1),d3(1),x0(:,1));
    
    if obj.iter==1 && strcmp(obj.modality,'MRI')
        % Change DC component of bias field to make intensities more
        % simillar between MR images. The scaling parameter is set in
        % init_buf
        b1 = chan(n).B1(1,1);
        b2 = chan(n).B2(1,1);
        b3 = chan(n).B3(1,1);
                
        chan(n).T(1,1,1) = 1/(b1*b2*b3)*log(scl(n));
    end
end

% Create initial bias field
%-----------------------------------------------------------------------
llrb = 0;
for n=1:N
    B1 = chan(n).B1;
    B2 = chan(n).B2;
    B3 = chan(n).B3;
    C  = chan(n).C;
    T  = chan(n).T;
    chan(n).ll = double(-0.5*T(:)'*C*T(:));
    for z=1:numel(z0)
        bf           = transf(B1,B2,B3(z,:),T);
        tmp          = bf(buf(z).msk{n}); 
        buf(z).bf{n} = single(exp(tmp));
    end
    llrb = llrb + chan(n).ll;
end

varargout{1} = buf;
varargout{2} = chan;
varargout{3} = llrb;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
return;
%=======================================================================