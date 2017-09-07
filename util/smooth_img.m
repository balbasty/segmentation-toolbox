function smooth_img(img,fwhm)
% Convolve the volume in memory (fwhm in voxels).
% FORMAT smooth_img(img,fwhm)
%__________________________________________________________________________

if nargin < 2, fwhm = [0.75 0.75 0.75]; end

lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(img,img,x,y,z,-[i j k]);
%==========================================================================