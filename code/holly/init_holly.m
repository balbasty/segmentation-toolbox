function holly = init_holly(pth_obj,dir_data,jnam)
HOLLY_CORES = 192;
RAM         = 6; % This value will be decreased automatically to an optimal value after the first iteration

holly = struct;

if numel(jnam)>9
    error('Name of Holly job needs to be less than 10 characters!')
end

holly.jnam_h     = jnam; 
holly.jnam_dummy = ['w' jnam];

[holly.dir_root_h,holly.dir_root_l,holly.dir_pwd_h] = read_directory_details('directory_details.txt',holly.jnam_h);    

holly.dir_matlab_h              = '/share/apps/MATLAB/R2016a/bin/matlab';
[holly.username,holly.password] = read_user_details('user_details.txt');   

holly.RAM = RAM;

M       = numel(pth_obj);
holly.t = 0;
for m=1:M
    S       = numel(pth_obj{m});
    holly.t = holly.t + S; 
end

split = ceil(holly.t/HOLLY_CORES);
fprintf('Splitting Holly array jobs into %d sub-job(s)\n\n',split);

holly.t = ceil(holly.t/split);

holly.dir_logs_l     = fullfile(holly.dir_root_l,'logs');  
holly.dir_scripts_l  = fullfile(holly.dir_root_l,'scripts');

if exist(holly.dir_logs_l,'dir'),     rmdir(holly.dir_logs_l,'s');     end; mkdir(holly.dir_logs_l);
if exist(holly.dir_scripts_l,'dir'),  rmdir(holly.dir_scripts_l,'s');  end; mkdir(holly.dir_scripts_l);

holly.dir_logs_h     = fullfile(holly.dir_root_h,'logs');           
holly.dir_scripts_h  = fullfile(holly.dir_root_h,'scripts');
      
holly.fun = 'segment_holly'; 

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

dir_holly_mat = fullfile(dir_data,'mat');
if exist(dir_holly_mat,'dir'), rmdir(dir_holly_mat,'s'); end; mkdir(dir_holly_mat);

pth_pth_obj = fullfile(dir_holly_mat,'pth_obj.mat');
save(pth_pth_obj,'pth_obj1') ;

pth_split = fullfile(dir_holly_mat,'split.mat');
save(pth_split,'split') ;
%==========================================================================