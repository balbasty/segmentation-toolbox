function [dir_root_h,dir_root_l,dir_pwd_h] = read_directory_details(fname,jnam)
fid = fopen(fname);

dir_root_h = fgetl(fid);
dir_root_l = fgetl(fid);
dir_pwd_h  = fgetl(fid);

dir_root_h = fullfile(dir_root_h,jnam);
dir_root_l = fullfile(dir_root_l,jnam);

fclose(fid);