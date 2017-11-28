function [holly,pth_pth_obj] = init_holly(pth_obj,dir_data,run_on_holly,holly)
jnam       = holly.jnam; 
jnam_dummy = holly.jnam_dummy;
RAM        = holly.RAM;
split      = holly.split;

holly       = struct;
pth_pth_obj = '';
if run_on_holly
    [holly.dir_root_h,holly.dir_root_l,holly.dir_pwd_h] = read_directory_details('directory_details.txt');
    holly.dir_matlab_h                                  = '/share/apps/MATLAB/R2016a/bin/matlab';
    [holly.username,holly.password]                     = read_user_details('user_details.txt');   

    holly.RAM = RAM;

    M       = numel(pth_obj);
    holly.t = 0;
    for m=1:M
        S       = numel(pth_obj{m});
        holly.t = holly.t + S; 
    end
    
    if split>holly.t
        split = min(holly.t,4);
    end
    
    holly.t = ceil(holly.t/split);
    
    holly.dir_logs_l     = fullfile(holly.dir_root_l,'logs');  
    holly.dir_scripts_l  = fullfile(holly.dir_root_l,'scripts');

    if exist(holly.dir_logs_l,'dir'),     rmdir(holly.dir_logs_l,'s');     end; mkdir(holly.dir_logs_l);
    if exist(holly.dir_scripts_l,'dir'),  rmdir(holly.dir_scripts_l,'s');  end; mkdir(holly.dir_scripts_l);
          
    holly.dir_logs_h     = fullfile(holly.dir_root_h,'logs');           
    holly.dir_scripts_h  = fullfile(holly.dir_root_h,'scripts');

    holly.jnam_h     = jnam;        
    holly.fun        = 'segment_holly'; 
    holly.jnam_dummy = jnam_dummy;      

    holly = create_bash_scripts(holly);   

    cnt      = 1;
    pth_obj1 = cell(1,holly.t);
    for m=1:M    
        S = numel(pth_obj{m});
        for s=1:S
            fname         = pth_obj{m}{s};
            nfname        = ['/' fname(7:end)];
            pth_obj1{cnt} = nfname;            
            cnt           = cnt + 1;
        end
    end

    pth_pth_obj = fullfile(dir_data,'pth_obj.mat');
    save(pth_pth_obj,'pth_obj1') ;
    
    pth_split = fullfile(dir_data,'split.mat');
    save(pth_split,'split') ;
end
%==========================================================================