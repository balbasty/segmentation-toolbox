function show_resp(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        K = numel(obj{m}{s}.pth_resp);
        
        for k=1:K    
            V  = spm_vol(obj{m}{s}.pth_resp{k});
            dm = V.dim;                        
            if dm(3)>1
                zix = floor(dm(3)/2) + 1;
                img = V.private.dat(:,:,zix);            
            else
                img = V.private.dat(:,:,1);            
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