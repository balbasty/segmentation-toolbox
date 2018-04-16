function show_resp(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        K = numel(obj{m}{s}.pth_resp);
        
        for k=1:K    
            Nii = nifti(obj{m}{s}.pth_resp{k});
            img = Nii.dat(:,:,:);
            dm  = size(img);
            if numel(dm)==3
                zix = floor(dm(3)/2) + 1;
                img = img(:,:,zix);            
            end
            
            wp = round(obj{m}{s}.segment.wp(k),3);
            
            subplot(M*numel(rand_subjs{1}),K,cnt);
            imagesc(img'); axis off xy; colormap(gray);
            title(['q_{' num2str(m), ',' num2str(s) ',' num2str(k) '},w=' num2str(wp)]);
            cnt = cnt + 1;
        end 
    end
end                                  
drawnow
%==========================================================================