function show_template(fig,pth_template)
% Display soft-maxed template
% FORMAT show_template(fig,pth_template)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

set(0,'CurrentFigure',fig);     

Nii = nifti(pth_template);
b   = spm_matcomp('softmax',Nii.dat(:,:,:,:));
d   = size(b);
Kb  = d(4);    
if d(3)>1
    for k=1:Kb         
        subplot(3,Kb,k);
        slice = b(:,:,floor(d(3)/2) + 1,k);
        imagesc(slice'); axis off xy; title(['k=' num2str(k)]); colormap(gray);               

        subplot(3,Kb,Kb + k);
        slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
        imagesc(slice); axis off xy; colormap(gray);   

        subplot(3,Kb,2*Kb + k);
        slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
        imagesc(slice'); axis off xy; colormap(gray);   
    end 
else
    K1 = floor(sqrt(Kb));
    K2 = ceil(Kb/K1);      
    for k=1:Kb
        subplot(K1,K2,k);
        slice = b(:,:,floor(d(3)/2) + 1,k);
        imagesc(slice'); axis xy off; title(['k=' num2str(k)]); colormap(gray);  
    end      
end
drawnow
%==========================================================================