function [V,labels] = get_V(im)

sanity_check_imobj(im);

S_requested = im{2};

[V,S,N] = read_ims_from_dir(im{1},S_requested);
labels  = read_ims_from_dir(im{5},S,N);

fprintf('Loaded data from %d subject(s) having %d channel(s) each\n',S,N); 
%==========================================================================

%==========================================================================
function [V,S,N] = read_ims_from_dir(pth,S,N)
if nargin<3, N=0; end

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

%==========================================================================
function sanity_check_imobj(im)
if numel(im)~=5
    error('numel(imobj)~=5')
end

folder    = dir(im{1}); 
folder    = folder(3:end);
dirflag   = [folder.isdir];

subfolder = folder(dirflag);
S         = numel(subfolder);
if S==0
    error('S==0')
end

files = dir(fullfile(im{1},subfolder(1).name,'*.nii'));
N0    = numel(files);
for s=2:S
    files = dir(fullfile(im{1},subfolder(s).name,'*.nii'));
    N1    = numel(files);
    if N0~=N1
        error('N0~=N1')
    end
end
%==========================================================================