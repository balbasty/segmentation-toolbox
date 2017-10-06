function [V,N,S] = get_V(imdir,S)
folder    = dir(imdir);         % folder with subfolders containing multi-channel data of subjects
folder    = folder(3:end);
dirflag   = [folder.isdir];
if sum(dirflag)==0
    % Check if folder contains Niftis
    N       = 1;
    files   = dir(fullfile(imdir,'*.nii'));
    S1      = numel(files);
    if S>S1
        S = S1;
    end
    if S
        V = cell(S,1);
        for s=1:S
            V{s} = spm_vol(fullfile(imdir,files(s).name));
        end
    else
        error('No Niftis in folder!')
    end
else
    subfolder = folder(dirflag);    % subfolders (S1,S2,...)
    S1        = numel(subfolder);
    if S>S1
        S = S1;
    end

    files = dir(fullfile(imdir,subfolder(1).name,'*.nii'));
    N     = numel(files);

    V = cell(S,1);
    for s=1:S
        folder = fullfile(imdir,subfolder(s).name);
        files  = dir(fullfile(folder,'*.nii'));
        for d=1:N
            V{s}(d) = spm_vol(fullfile(folder,files(d).name));
        end
    end
end

fprintf('Loading data from %d subjects having %d channels each\n',S,N); 
%==========================================================================