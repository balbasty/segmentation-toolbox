function slice_template(fname)
if nargin<2, axis_2d = 3; end

V  = spm_vol(fname);
dm = V(1).dim;
K  = numel(V);

if axis_2d==1
    d1 = floor(dm(1)/2) + 1;
    bb = [d1 d1;-inf inf;-inf inf];
elseif axis_2d==2
    d1 = floor(dm(2)/2) + 1;
    bb = [-inf inf;d1 d1;-inf inf];
elseif axis_2d==3 
    d1 = floor(dm(3)/2) + 1;
    bb = [-inf inf;-inf inf;d1 d1];
end                

for k=1:K
    spm_impreproc('subvol',V(k),bb','2d_');      
end