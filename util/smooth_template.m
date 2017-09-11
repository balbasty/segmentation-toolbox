function [munum,muden] = smooth_template(munum,muden,d,fwhm)
if nargin<4, fwhm = 0.5; end

if ~fwhm
    return;
end

fwhm = fwhm*ones(1,3);

lim = ceil(2*fwhm);
x   = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x = x/sum(x);
y   = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y = y/sum(y);
z   = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z = z/sum(z);
i   = (length(x) - 1)/2;
j   = (length(y) - 1)/2;
k   = (length(z) - 1)/2;

K     = size(munum,1);
munum = reshape(munum',[d K]);
muden = reshape(muden',[d K]);
for k1=1:K
    img = munum(:,:,:,k1);            
    spm_conv_vol(img,img,x,y,z,-[i j k]);
    munum(:,:,:,k1) = img;

    img = muden(:,:,:,k1);            
    spm_conv_vol(img,img,x,y,z,-[i j k]);
    muden(:,:,:,k1) = img;
end
munum = reshape(munum,[prod(d) K])';
muden = reshape(muden,[prod(d) K])';