function [V,labels] = get_V(im)

sanity_check_imobj(im);

S_requested = im{2};

[V,S,N] = read_ims_from_dir(im{1},S_requested);
labels  = read_ims_from_dir(im{6},S,N);

fprintf('Loaded data from %d subject(s) having %d channel(s) each\n',S,N); 
%==========================================================================

%==========================================================================
function sanity_check_imobj(im)
if numel(im)~=6
    error('numel(imobj)~=6')
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
        subfolder(s).name
        error('N0~=N1')
    end
end
%==========================================================================