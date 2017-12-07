function RAM = estimate_RAM_usage(jobid,oRAM)
if nargin<2, oRAM = 6; end

[holly.username,holly.password] = read_user_details('user_details.txt');   

cmd             = ['sshpass -p "' holly.password '" ssh -o StrictHostKeyChecking=no ' holly.username '@holly "source /etc/profile;/opt/gridengine/bin/linux-x64/qacct -j ' num2str(jobid) ' | grep maxvmem"'];
[status,result] = system(cmd);   

if status==0
    jobs = strsplit(result,'G');
    S    = numel(jobs) - 1;    
    a    = zeros(1,S);
    for s=1:S
        job  = jobs{s};                                     % string
        val  = str2double(regexp(job,'\d','match','once')); % find value of first digit in string
        ix   = strfind(job,num2str(val));                   % index of first digit in string
        a(s) = str2double(job(ix(1):end));                  % extract RAM
    end
    
    sd  = 0.2;
    RAM = round(max(a) + sd*max(a),1);

    fprintf('New RAM usage on holly: %d GB\n',RAM);
else
    fprintf([result '\n'])
    RAM = oRAM;
    warning('status~=0, RAM set to previous RAM')    
end
