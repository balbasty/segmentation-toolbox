function [s,ds1,ds2,ds3] = spm_sample_logpriors(tpm,x1,x2,x3)
deg = tpm.deg;
wp  = tpm.wp;
d   = tpm.V.dim;
dx  = size(x1);
Kb  = numel(tpm.dat);

% For masking outside FOV
msk1 = x1>=1 & x1<=d(1) & x2>=1 & x2<=d(2) & x3>=1 & x3<=d(3);
msk2 = x3<1;
x1   = x1(msk1);
x2   = x2(msk1);
x3   = x3(msk1);

s = cell(1,Kb);
if nargout<=1    
    mx = -Inf(dx);
    for k=1:Kb
        a = spm_bsplins(tpm.dat{k},x1,x2,x3,[deg deg deg 0 0 0]); % interpolate log-template
                        
        s{k}       = ones(dx)*log(wp(k));
        s{k}(msk1) = a;
        s{k}(msk2) = log(wp(k));     
        
        mx = max(mx,s{k});
    end
    
    % Safe soft-max
    tot = zeros(dx);
    for k=1:Kb
        s{k} = exp(s{k} - mx);
        tot  = tot + s{k}; 
    end    
    
    for k=1:Kb
        s{k} = s{k}./tot;
    end
else
    mx  = -Inf(dx);
    ds1 = cell(1,Kb);
    ds2 = cell(1,Kb);
    ds3 = cell(1,Kb);    
    for k=1:Kb
        [a,da1,da2,da3] = spm_bsplins(tpm.dat{k},x1,x2,x3,[deg deg deg 0 0 0]); % interpolate log-template and compute derivatives     
                
        s{k}       = ones(dx)*log(wp(k));
        s{k}(msk1) = a;
        s{k}(msk2) = log(wp(k));        
        
        mx = max(mx,s{k});
        
        ds1{k} = zeros(dx); ds1{k}(msk1) = da1;
        ds2{k} = zeros(dx); ds2{k}(msk1) = da2;
        ds3{k} = zeros(dx); ds3{k}(msk1) = da3;
    end
    
    % Safe soft-max
    tot = zeros(dx);
    for k=1:Kb
        s{k} = exp(s{k} - mx);
        tot  = tot + s{k}; 
    end    
    
    da1 = zeros(dx);
    da2 = zeros(dx);
    da3 = zeros(dx);
    for k=1:Kb
         s{k} = s{k}./tot;
         da1  = da1 + s{k}.*ds1{k};
         da2  = da2 + s{k}.*ds2{k};
         da3  = da3 + s{k}.*ds3{k};
    end
    for k=1:Kb
        ds1{k} = s{k}.*(ds1{k} - da1);
        ds2{k} = s{k}.*(ds2{k} - da2);
        ds3{k} = s{k}.*(ds3{k} - da3);
    end
end