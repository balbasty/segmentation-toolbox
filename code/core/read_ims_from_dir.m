function [V,S,N] = read_ims_from_dir(pth,S,N)
if nargin<2, S = Inf; end
if nargin<3, N = 1; end

if isempty(pth)
    % No labels provided
    V = cell(1,S);
    for s=1:S
        for n=1:N
            V{s}(n) = NaN;
        end
    end
else
    folder    = dir(pth);  % folder with subfolders containing multi-channel data of subjects
    folder    = folder(3:end);
    dirflag   = [folder.isdir];

    subfolder = folder(dirflag);   % subfolders (S1,S2,...)
    S1        = numel(subfolder);
    if S>S1
        S = S1;
    end

    files = dir(fullfile(pth,subfolder(1).name,'*.nii'));
    N     = numel(files);

    V = cell(1,S);
    for s=1:S
        folder = fullfile(pth,subfolder(s).name);
        files  = dir(fullfile(folder,'*.nii'));
        for n=1:N
            V{s}(n) = spm_vol(fullfile(folder,files(n).name));
        end        
    end 
end
%==========================================================================