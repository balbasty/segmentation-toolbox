function [V,N,totS] = get_V(imdir)
totS = 0;
for i=1:numel(imdir)
    S_requested = imdir{i}{3};
    
    folder    = dir(imdir{i}{1});         % folder with subfolders containing multi-channel data of subjects
    folder    = folder(3:end);
    dirflag   = [folder.isdir];

    subfolder = folder(dirflag);    % subfolders (S1,S2,...)
    S1        = numel(subfolder);
    if S_requested>S1
        S = S1;
    else
        S = S_requested;
    end

    files = dir(fullfile(imdir{i}{1},subfolder(1).name,'*.nii'));
    N     = numel(files);

    for s=1:S
        totS   = totS + 1;
        folder = fullfile(imdir{i}{1},subfolder(s).name);
        files  = dir(fullfile(folder,'*.nii'));
        for d=1:N
            V{totS}(d) = spm_vol(fullfile(folder,files(d).name));
        end        
    end 
end
fprintf('Loading data from %d subjects having %d channels each\n',totS,N); 
%==========================================================================