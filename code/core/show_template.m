function show_template(pth_template,fig)
if nargin<2, fig = figure(666); end

set(0,'CurrentFigure',fig);     

Nii = nifti(pth_template);
b   = exp(Nii.dat(:,:,:,:));
b   = bsxfun(@rdivide,b,sum(b,4));
d   = size(b);
Kb  = d(4);    
if d(3)>1
    for k=1:Kb         
        subplot(3,Kb,k);
        slice = b(:,:,floor(d(3)/2) + 1,k);
        imagesc(slice'); axis off image xy; title(['k=' num2str(k)]); colormap(gray);               

        subplot(3,Kb,Kb + k);
        slice = permute(b(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
        imagesc(slice); axis off image xy; colormap(gray);   

        subplot(3,Kb,2*Kb + k);
        slice = permute(b(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
        imagesc(slice'); axis off image xy; colormap(gray);   
    end 
else
    K1 = floor(sqrt(Kb));
    K2 = ceil(Kb/K1);      
    for k=1:Kb
        subplot(K1,K2,k);
        slice = b(:,:,floor(d(3)/2) + 1,k);
        imagesc(slice'); axis image xy off; title(['k=' num2str(k)]); colormap(gray);  
    end      
end
drawnow
%==========================================================================