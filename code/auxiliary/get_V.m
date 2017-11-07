function [V,S,N] = get_V(imobj)

sanity_check_imobj(imobj);

S_requested = imobj{2};

folder    = dir(imobj{1});  % folder with subfolders containing multi-channel data of subjects
folder    = folder(3:end);
dirflag   = [folder.isdir];

subfolder = folder(dirflag);   % subfolders (S1,S2,...)
S1        = numel(subfolder);
if S_requested>S1
    S = S1;
else
    S = S_requested;
end

files = dir(fullfile(imobj{1},subfolder(1).name,'*.nii'));
N     = numel(files);

V = cell(1,S);
for s=1:S
    folder = fullfile(imobj{1},subfolder(s).name);
    files  = dir(fullfile(folder,'*.nii'));
    for n=1:N
        V{s}(n) = spm_vol(fullfile(folder,files(n).name));
    end        
end 
fprintf('Loaded data from %d subjects having %d channels each\n',S,N); 
%==========================================================================

%==========================================================================
function sanity_check_imobj(imobj)
if numel(imobj)~=4
    error('numel(imobj)~=4')
end

folder    = dir(imobj{1}); 
folder    = folder(3:end);
dirflag   = [folder.isdir];

subfolder = folder(dirflag);
S         = numel(subfolder);
if S==0
    error('S==0')
end

files = dir(fullfile(imobj{1},subfolder(1).name,'*.nii'));
N0    = numel(files);
for s=2:S
    files = dir(fullfile(imobj{1},subfolder(s).name,'*.nii'));
    N1    = numel(files);
    if N0~=N1
        error('N0~=N1')
    end
end
%==========================================================================