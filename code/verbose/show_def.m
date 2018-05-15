function show_def(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        V   = spm_vol(obj{m}{s}.pth_vel);
        dm  = V.dim;
        zix = floor(dm(3)/2) + 1;                                         
        
        subplot(M,numel(rand_subjs{1}),cnt)
        cnt = cnt + 1;
        
        img = [];
        for i=1:3                        
            slice = single(V(1).private.dat(:,:,zix,i));            
            img   = [img slice'];            
        end 
        
        imagesc(img); axis off image xy; colormap(gray);
    end
end                                  
drawnow
%==========================================================================