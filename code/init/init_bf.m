function [buf,chan,llrb] = init_bf(buf,N,obj,V,x0,y0,z0,ff)
kron      = @(a,b) spm_krutil(a,b);
avg_bf_dc = obj.avg_bf_dc;
tot_S     = obj.tot_S;
ml        = obj.do_ml;

cl   = cell(N,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl};
chan = struct(args{:});
for n=1:N
    biasreg = obj.biasreg(n)*ff; 
    vx      = vxsize(V(n).mat);
    fwhm    = obj.biasfwhm(n);
    d0      = V(n).dim;
    
    sd = vx(1)*d0(1)/fwhm; d3(1) = ceil(sd*2); krn_x = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
    sd = vx(2)*d0(2)/fwhm; d3(2) = ceil(sd*2); krn_y = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
    sd = vx(3)*d0(3)/fwhm; d3(3) = ceil(sd*2); krn_z = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
    
    if tot_S==1
        % GAUSSIAN REGULARISATION for bias correction
        Cbias = kron(krn_z,kron(krn_y,krn_x)).^(-2)*biasreg;        
        chan(n).C = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias)); % Store prior covaricance for bias regularisation
    else
        % BENDING ENERGY regularisation for bias correction
        % This penalises the sum of squares of the 2nd derivatives of the bias parameters
        chan(n).C = diffeo('penalty',[numel(krn_x) numel(krn_y) numel(krn_z)],vx,[0 0 biasreg 0 0]);
        chan(n).C = chan(n).C(1:size(chan(n).C,1)/3,1:size(chan(n).C,1)/3);
    end
    
    % Initial parameterisation of bias field
    if isfield(obj,'Tbias') && ~isempty(obj.Tbias{n})                
        chan(n).T                    = obj.Tbias{n};
        if tot_S>1, chan(n).T(1,1,1) = chan(n).T(1,1,1) - avg_bf_dc(n); end
    else
        chan(n).T = zeros(d3);
    end
    
    % Basis functions for bias correction
    chan(n).B3 = spm_dctmtx(d0(3),d3(3),z0);
    chan(n).B2 = spm_dctmtx(d0(2),d3(2),y0(1,:)');
    chan(n).B1 = spm_dctmtx(d0(1),d3(1),x0(:,1));
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
        bf                  = transf(B1,B2,B3(z,:),T);
        tmp                 = bf(buf(z).msk{n});
        if ml, chan(n).ll   = chan(n).ll + double(sum(tmp)); end        
        buf(z).bf{n}        = single(exp(tmp));
    end
    llrb = llrb + chan(n).ll;
    clear B1 B2 B3 T C
end
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