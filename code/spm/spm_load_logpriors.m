function tpm = spm_load_logpriors(pth,wp)
V = spm_vol(pth); 
spm_check_orientations(V);

tpm.deg = 1; % Degree of interpolation
tpm.V   = V;
tpm.M   = tpm.V(1).mat;
tpm.wp  = wp;

% Allocate data
Kb = numel(tpm.V);
tpm.dat = cell(Kb,1);
for k1=1:(Kb)
    tpm.dat{k1} = zeros(tpm.V(1).dim(1:3));
end

% Slice log-template
for i=1:tpm.V(1).dim(3)
    M = spm_matrix([0 0 i]);
    for k1=1:Kb
        tpm.dat{k1}(:,:,i) = spm_slice_vol(tpm.V(k1),M,tpm.V(1).dim(1:2),0);
    end
end

% Convert into b-spline coefficients
tpm.bg1 = zeros(Kb,1);
tpm.bg2 = zeros(Kb,1);
for k1=1:Kb
    tpm.bg1(k1) = mean(mean(tpm.dat{k1}(:,:,1)));
    tpm.bg2(k1) = mean(mean(tpm.dat{k1}(:,:,end)));
    tpm.dat{k1} = spm_bsplinc(tpm.dat{k1},[tpm.deg tpm.deg tpm.deg 0 0 0]);
end