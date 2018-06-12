function obj = modify_obj(obj,iter,pars)
sparam0 = pars.sparam0;
sparam1 = pars.sparam1;
it1     = numel(fliplr(sparam1:10:sparam0)) + 1;
niter   = pars.niter;
M       = numel(obj);    
for m=1:M
    S = numel(obj{m});    
    for s=1:S            
        obj{m}{s}.iter = iter;
        
        if iter==2            
            obj{m}{s}.uniform = false;  
            
            obj{m}{s}.segment.nitgmm = 10;       
            
            % Enable updating bias field
            obj{m}{s}.segment.do_bf  = obj{m}{s}.segment.do_bf0;                                                                                                 
        end

        if iter>=2
            % Decrease registration regularisation over iterations
            reg0  = obj{m}{s}.segment.reg0;   
            sched = 2.^fliplr(1:8);
            
            if iter>=it1                
                scl = sched(min((iter + 1) - it1,numel(sched)));               
            else
                scl = sched(1);
            end
            
            obj{m}{s}.segment.reg(3) = reg0(3)*scl;       
        end        
        
        if iter==it1
            % Enable updating tissue and Gaussian weights
            obj{m}{s}.segment.do_wp = obj{m}{s}.segment.do_wp0;   
            obj{m}{s}.segment.do_mg = obj{m}{s}.segment.do_mg0;   
            
            % Enable updating MRF weights
            obj{m}{s}.segment.mrf.update_Upsilon = obj{m}{s}.segment.mrf.update_Upsilon0;
        end                
        
        if iter==niter
            % Final iteration -> write some output
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