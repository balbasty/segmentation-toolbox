function show_template(fig,pth_template,dir_figs,iter)
% Display soft-maxed template
% FORMAT show_template(fig,pth_template,dir_figs,iter)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

set(0,'CurrentFigure',fig); clf(fig);

Nii = nifti(pth_template);
b   = spm_matcomp('softmax',Nii.dat(:,:,:,:));
d   = size(b);
Kb  = d(4);    
if d(3)>1
    img1 = [];
    img2 = [];
    img3 = [];
    for k=1:Kb         
        slice = b(:,:,floor(d(3)/2) + 1,k);
        img1  = [img1 slice'];
        
        slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
        img2  = [img2 slice];
        
        slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
        img3  = [img3 slice'];
    end 
    
    subplot(311);
    imagesc(img1,[0 1]); axis off image xy; colormap(gray);      
    title(['Template (iter=' num2str(iter) ')'])

    subplot(312);
    imagesc(img2,[0 1]); axis off image xy; colormap(gray);   

    subplot(313);
    imagesc(img3,[0 1]); axis off image xy; colormap(gray);   
else    
    img = [];
    for k=1:Kb        
        slice = b(:,:,floor(d(3)/2) + 1,k);
        img   = [img slice'];        
    end      
    
    imagesc(img,[0 1]); axis xy off image; colormap(gray);     
    title(['template (iter=' num2str(iter) ')'])
end
drawnow
%==========================================================================