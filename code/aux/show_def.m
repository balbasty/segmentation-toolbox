function show_def(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        Nii = nifti(obj{m}{s}.pth_vel);
        img = Nii.dat(:,:,:,:);
        dm  = size(img);
        zix = floor(dm(3)/2) + 1; 
        
        for i=1:3
            subplot(M*numel(rand_subjs{1}),3,cnt)
            imagesc(img(:,:,zix,i)'); axis off xy; colormap(gray); colorbar
            title(['def_{' num2str(m), ',' num2str(s), ',' num2str(i) '}']);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================