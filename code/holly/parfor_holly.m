function holly = parfor_holly(holly)

% Submit subject specific jobs (i.e., estimate on all subjects)
pause(1);
tic;
cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qsub -l vf=' num2str(holly.RAM) 'G -l h_vmem=' num2str(holly.RAM) 'G ' holly.pth_script_parfor '"'];        
[status,result] = system(cmd);    
if status
    fprintf([result '\n'])
    error('status~=0 on Holly!') 
end
fprintf(result)

jobid = result(15:20); % ID of array job on Holly

% Submit dummy job
pause(1);
cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qsub -l vf=0.1G -l h_vmem=0.1G -hold_jid ' holly.jnam_h ' -cwd ' holly.pth_script_dummy '"'];
[status,result] = system(cmd);
if status
    fprintf([result '\n'])
    error('status~=0 for dummy job on Holly!') 
end
fprintf(result)

cmd = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qstat | grep ' holly.jnam_dummy '"'];        
while 1, 
    % Check if dummy job has finished
    pause(2);

    [status,result] = system(cmd);   
    if isempty(result)
        fprintf('Elapsed time (holly): %d s\n',round(toc))                        
        pause(1);
        break
    end
end  

holly.RAM = estimate_RAM_usage(jobid,holly.RAM);
%========================================================================== 