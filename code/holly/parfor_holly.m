function [ll,munum,muden,Nm] = parfor_holly(pth_obj,holly,num_workers)
% Submit subject specific jobs
tic;
cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qsub -l vf=' num2str(holly.RAM) 'G -l h_vmem=' num2str(holly.RAM) 'G ' holly.pth_script_parfor '"'];        
[status,result] = system(cmd);    
if status
    fprintf([result '\n'])
    error('status~=0 on Holly!') 
end
fprintf(result)

% Submit dummy job
cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qsub -l vf=0.1G -l h_vmem=0.1G -hold_jid ' holly.jnam_h ' -cwd ' holly.pth_script_dummy '"'];
[status,result] = system(cmd);
if status
    fprintf([result '\n'])
    error('status~=0 for dummy job on Holly!') 
end
fprintf(result)

while 1, 
    % Check if dummy job has finished
    pause(1);

    cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qstat | grep ' holly.jnam_dummy '"'];        
    [status,result] = system(cmd);   
    if isempty(result)
        fprintf('Elapsed time (holly): %d s\n',round(toc))
                
        M     = numel(pth_obj);
        munum = 0; muden = 0; ll = 0; Nm = 0;
        for m=1:M
            S         = numel(pth_obj{m});
            pth_obj_m = pth_obj{m};                
%                 for s=1:S
            manage_parpool(num_workers);
            parfor (s=1:S,num_workers)
                obj = load(pth_obj_m{s},'-mat','munum','muden','ll','nm');

                munum = munum + obj.munum;
                muden = muden + obj.muden;
                ll    = ll    + obj.ll;
                Nm    = Nm    + obj.nm;
            end
        end

        return
    end
end  
%========================================================================== 