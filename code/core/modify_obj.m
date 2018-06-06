function obj = modify_obj(obj,iter,niter)
M = numel(obj);    
for m=1:M
    S = numel(obj{m});    
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==2            
            obj{m}{s}.uniform = false;  
            
            obj{m}{s}.segment.do_bf  = obj{m}{s}.segment.do_bf0;                                                                                     
            obj{m}{s}.segment.niter  = 3;
            obj{m}{s}.segment.nitgmm = 20;                                        
        end

        if iter==3
            obj{m}{s}.segment.mrf.update_Upsilon = obj{m}{s}.segment.mrf.update_Upsilon0;
        end
        
        if iter==6
            obj{m}{s}.segment.do_wp = obj{m}{s}.segment.do_wp0;   
            obj{m}{s}.segment.do_mg = obj{m}{s}.segment.do_mg0;   
        end
        
        if iter>=2
            reg0                     = obj{m}{s}.segment.reg0;   
            sched                    = 2.^fliplr(repelem(0:9,2));
            scl                      = sched(min(iter,numel(sched)));               
            obj{m}{s}.segment.reg(3) = reg0(3)*scl;       
        end
        
        if iter==niter
            obj{m}{s}.write_res.do_write_res  = true; 
            obj{m}{s}.write_res.write_tc(:,:) = true;             
            obj{m}{s}.write_res.write_bf(:,:) = true;             
            obj{m}{s}.write_res.write_df(:,:) = true;       
            
            if obj{m}{s}.image(1).dim(3)==1
                % Is 2D
                obj{m}{s}.write_res.write_tc(:,2) = false; 
                obj{m}{s}.write_res.write_df(:,:) = false; 
            end
        end
    end
end
%==========================================================================