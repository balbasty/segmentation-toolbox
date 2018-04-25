function show_subj_ll(fig,obj,rand_subjs)
set(0,'CurrentFigure',fig);       

M   = numel(obj);
cnt = 1;    
for m=1:M                
    for s=rand_subjs{m} 
        subplot(M,numel(rand_subjs{1}),cnt)
        
        ll = obj{m}{s}.segment.ll_all;        
                
        plot(0:numel(ll) - 1,ll,'b-'); title(['ll_{m=' num2str(m) ', s=' num2str(s) '}'])
        
        cnt = cnt + 1;
    end
end                                  
drawnow
%==========================================================================