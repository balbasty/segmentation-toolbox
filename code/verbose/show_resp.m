function show_resp(fig,obj,rand_subjs,iter,print_wp)
if nargin<5, print_wp = false; end

set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        K = numel(obj{m}{s}.pth_resp);
        
        subplot(M*numel(rand_subjs{1}),1,cnt);
        
        nam   = 'wp = [';
        img = [];
        for k=1:K    
            V  = spm_vol(obj{m}{s}.pth_resp{k});
            dm = V.dim;                        
            if dm(3)>1
                zix = floor(dm(3)/2) + 1;
                slice = V.private.dat(:,:,zix);            
            else
                slice = V.private.dat(:,:,1);            
            end
            img = [img slice'];
            
            wp = round(obj{m}{s}.segment.wp(k),2);  
            if k==1
                nam = [nam num2str(wp)];
            elseif k==K
                nam = [nam ', ' num2str(wp) ']'];                
            else
                nam = [nam ', ' num2str(wp)];
            end            
        end 
        
        imagesc(img); axis off xy; colormap(gray);
        if cnt==1 && ~print_wp
            title(['Template-space responsibilities (iter=' num2str(iter) ')'])           
        elseif print_wp
            title(nam);            
        end
        
        cnt = cnt + 1;
    end
end                                  
drawnow
%==========================================================================