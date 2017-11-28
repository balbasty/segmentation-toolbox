function holly = create_bash_scripts(holly)
% Create bash scripts that are called using qsub to run on holly.
%
% Mikael Brudfors
% 2017-11-04
%==========================================================================

dir_scripts_l = holly.dir_scripts_l;
dir_scripts_h = holly.dir_scripts_h;
jnam_h        = holly.jnam_h;
jnam_dummy    = holly.jnam_dummy;
fun           = holly.fun;
dir_logs_h    = holly.dir_logs_h;
matlab_h      = holly.dir_matlab_h;
pwd_h         = holly.dir_pwd_h;
t             = holly.t;

% Create a bash-script defining what task should be parallelised on holly
bash_script_parfor = sprintf(['#!/bin/sh\n'...
                              '\n'...
                              '#$ -S /bin/sh\n'...
                              '#$ -N ' jnam_h '\n'...
                              '#$ -o ' fullfile(dir_logs_h,[jnam_h '_o.txt']) '\n'...
                              '#$ -e ' fullfile(dir_logs_h,[jnam_h '_e.txt']) '\n'...
                              '#$ -j n \n'...
                              '#$ -t 1-' num2str(t) ' \n'...
                              '\n'...
                              'matlab_command="addpath(genpath(''' pwd_h '''));' fun '($SGE_TASK_ID);quit"\n'...
                              matlab_h ' -nojvm -nodesktop -nosplash -nodisplay -singleCompThread -r $matlab_command']);

pth_script_parfor = fullfile(dir_scripts_l,[jnam_h '.sh']);
fid               = fopen(pth_script_parfor,'w');
fprintf(fid,bash_script_parfor);
fclose(fid);

fileattrib(pth_script_parfor,'+x','u')

pth_script_parfor = fullfile(dir_scripts_h,[jnam_h '.sh']);

% Create a bash-script defining a dummy task used a flag for detecting when
% the jobs running on holly have finished
bash_script_dummy = sprintf(['#!/bin/sh\n'...
                             '\n'...
                             '#$ -S /bin/sh\n'...
                             '#$ -N ' jnam_dummy '\n'...
                             '#$ -o ' fullfile(dir_logs_h,[jnam_dummy '_o.txt']) '\n'...
                             '#$ -e ' fullfile(dir_logs_h,[jnam_dummy '_e.txt']) '\n'...
                             '#$ -j n \n'...
                             '#$ -t 1-1 \n'...
                             '\n'...
                             'echo dummy\n']);

pth_script_dummy = fullfile(dir_scripts_l,[jnam_dummy '.sh']);
fid              = fopen(pth_script_dummy,'w');
fprintf(fid,bash_script_dummy);
fclose(fid);

fileattrib(pth_script_dummy,'+x','u')

pth_script_dummy = fullfile(dir_scripts_h,[jnam_dummy '.sh']);

holly.pth_script_parfor = pth_script_parfor;
holly.pth_script_dummy  = pth_script_dummy;
%==========================================================================