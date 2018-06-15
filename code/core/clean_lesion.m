function [P] = clean_lesion(P,level)
if nargin<2, level = 1; end

ix_bg  = 1; % 
ix_sk  = 2; %
ix_csf = 3; %
ix_sin = 4; % 
ix_st  = 5; %
ix_les = 6; %
ix_cal = 7; % 
ix_gm1 = 8; %
ix_gm2 = 9; %
ix_wm  = 10; % 

ix_gm = [ix_gm1 ix_gm2];

%--------------------------------------------------------------------------
b = sum(P(:,:,:,ix_les),4);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations of stuff inside of CSF
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp   = double(sum(P(:,:,i,ix_gm),4));
        wp   = double(sum(P(:,:,i,ix_wm),4));
        cp   = double(sum(P(:,:,i,ix_csf),4));
        lesp = double(sum(P(:,:,i,ix_les),4));
        
        bp = double(b(:,:,i))/255;
        bp = (bp>th).*(wp + gp + lesp + cp);
        
        b(:,:,i) = uint8(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end
clear gp wp lesp bp cp

%--------------------------------------------------------------------------
th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4)
        slices{k1} = double(P(:,:,i,k1))/255;
    end
    bp        = double(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{ix_gm1} + slices{ix_gm2} + slices{ix_wm} + slices{ix_sin} + slices{ix_les} + slices{ix_cal}))>th;
    slices{ix_gm1} = slices{ix_gm1}.*bp;
    slices{ix_gm2} = slices{ix_gm2}.*bp;
    slices{ix_wm}  = slices{ix_wm}.*bp;
    slices{ix_sin} = slices{ix_sin}.*bp;
    slices{ix_les} = slices{ix_les}.*bp;
    slices{ix_cal} = slices{ix_cal}.*bp;
    
    %----------------------------------------------------------------------
    if niter2>0
        cp             = double(c(:,:,i))/255;
        cp             = ((cp>th).*(slices{ix_csf} + slices{ix_gm1} + slices{ix_gm2} + slices{ix_wm} + slices{ix_sin} + slices{ix_les} + slices{ix_cal}))>th;
        slices{ix_csf} = slices{ix_csf}.*cp;
    end
    
    %----------------------------------------------------------------------
    tot = zeros(size(bp)) + eps;
    for k1=1:size(P,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4)
        P(:,:,i,k1) = uint8(round(slices{k1}./tot*255));
    end 
end
spm_progress_bar('Clear');
%==========================================================================