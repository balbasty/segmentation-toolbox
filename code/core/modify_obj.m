function obj = modify_obj(obj,iter,niter)
M = numel(obj);    
for m=1:M
    S         = numel(obj{m});    
    modality  = obj{m}{1}.modality;
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==1             
            obj{m}{s}.push_resp.do_push_resp = true;
            obj{m}{s}.write_res.do_write_res = false;
            
            obj{m}{s}.segment.do_def = false;
            obj{m}{s}.segment.do_bf  = false;                        
            obj{m}{s}.segment.niter  = 1;    
            
            obj{m}{s}.segment.do_wp  = false;
            obj{m}{s}.segment.nsubit = 1;
            obj{m}{s}.segment.nitgmm = 1;
        end

        if iter==2            
            obj{m}{s}.uniform = false;  
            
            obj{m}{s}.segment.do_bf  = obj{m}{s}.segment.do_bf0;                                                               
            
            obj{m}{s}.segment.do_wp  = obj{m}{s}.segment.do_wp0;   
            obj{m}{s}.segment.nsubit = 8;
            obj{m}{s}.segment.nitgmm = 20;  
        end

        if iter>=2
            reg0  = obj{m}{s}.segment.reg0;   
            sched = 2.^fliplr(repelem(0:9,2));
            scal  = sched(min(iter,numel(sched)));   
            
            obj{m}{s}.segment.reg(3) = reg0(3)*scal;                     
        end
        
        if iter==niter && obj{m}{s}.image(1).dim(3)>1
            obj{m}{s}.write_res.do_write_res = true;
            obj{m}{s}.write_res.mrf = 2;
            obj{m}{s}.write_res.write_tc(:,[1 2 4]) = true;            
        end
    end
end
%==========================================================================