function niis2subfolders(dir_niis,dir_subfolders)
files = dir(fullfile(dir_niis,'*.nii'));
S     = numel(files);

if exist(dir_subfolders,'dir'), rmdir(dir_subfolders,'s'); end; mkdir(dir_subfolders);

for s=1:S
    dir_s = fullfile(dir_subfolders,['S' num2str(s)]);
    mkdir(dir_s);
    
    fname = fullfile(dir_niis,files(s).name);
    
    copyfile(fname,dir_s);
end

