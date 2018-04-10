function obj = modify_obj(obj,iter,niter)
M = numel(obj);    
for m=1:M
    S         = numel(obj{m});    
    sum_bf_dc = 0;
    cnt_S     = 0;
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==1             
            obj{m}{s}.push_resp.do_push_resp = true;
            obj{m}{s}.write_res.do_write_res = false;
            
            obj{m}{s}.segment.do_def = false;
            obj{m}{s}.segment.do_bf  = false;
            obj{m}{s}.segment.do_wp  = false;
            obj{m}{s}.segment.niter  = 1;                 
            obj{m}{s}.segment.nsubit = 1;
            obj{m}{s}.segment.nitgmm = 1;
        end

        if iter==2            
            obj{m}{s}.uniform = false;  
            
            obj{m}{s}.segment.nsubit = 8;
            obj{m}{s}.segment.nitgmm = 20;  
            obj{m}{s}.segment.do_bf  = obj{m}{s}.segment.do_bf0;                                                  
            obj{m}{s}.segment.do_wp  = obj{m}{s}.segment.do_wp0;    
        end

        if iter>=2
            reg0  = obj{m}{s}.segment.reg0;   
            sched = 2.^fliplr(repelem(0:9,2));
            scal  = sched(min(iter,numel(sched)));   
            
            obj{m}{s}.segment.reg(3) = reg0(3)*scal;                     
        end
        
        % Sum bias field DC components
        sum_bf_dc = sum_bf_dc + obj{m}{s}.segment.bf_dc;
        cnt_S     = cnt_S + 1;
        
        if iter==niter && obj{m}{s}.image(1).dim(3)>1
            obj{m}{s}.write_res.do_write_res = true;
            obj{m}{s}.write_res.mrf = 2;
            obj{m}{s}.write_res.write_tc(:,[1 2 4]) = true;            
        end
    end
    
    % Average of bias field DC components
    avg_bf_dc = sum_bf_dc/cnt_S;
    
    % Set average bias field DC component
    for s=1:S 
        obj{m}{s}.segment.avg_bf_dc = avg_bf_dc; 
    end
end
%==========================================================================