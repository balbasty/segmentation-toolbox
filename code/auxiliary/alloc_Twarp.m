function pthTwarp = alloc_Twarp(obj,s)
d0 = obj.image(1).dim;
vx = sqrt(sum(obj.image(1).mat(1:3,1:3).^2));
sk = max([1 1 1],round(obj.samp*[1 1 1]./vx));
x0 = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0 = 1:sk(3):d0(3);
d  = [size(x0) length(z0)];

Twarp    = zeros([d 3],'single');    
pthTwarp = fullfile(obj.dirTwarp,['Twarp' num2str(s) '.mat']);     
save(pthTwarp,'Twarp');

