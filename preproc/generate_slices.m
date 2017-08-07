clear; close all; clc;

%% Get input data
N = 100;

datadir = '/home/mbrud/Data/IXI-clean/IXI-T1/';
pth     = dir(fullfile(datadir,'*.nii'));  

pthx = cell(N,1);
for n=1:N
    pthx{n,1} = fullfile(datadir,pth(n).name);
end

%% Align to MNI-space
addpath('../util/')
addpath('../core/')

fprintf('Preprocessing image: '); 
parfor n=1:N
    fprintf('%d ',n); 
     
    P{n} = ['im' num2str(n) '.nii'];
   
    cpy_img(pthx{n,1},P{n});                                  

    nm_reorient(P{n},1,1);
    reset_origin(P{n});

    rigid_align(P{n});  
end    
fprintf('\n'); 

%% Warp same size
samp = 1.5;
P    = warp_im_same_size(P,samp,'.');

%% Make 2D
prefix = 'slice-';

for n=1:N
    V  = spm_vol(P{n});
    d  = V.dim;
    d1 = floor(d(3)/2) + 25;
    subvol(V,[-inf inf;-inf inf;d1 d1]',prefix);
end

%% Visualise
NN = N;

F1 = floor(sqrt(NN));
F2 = ceil(NN/F1);      
for n=1:NN
    [dir,nam,ext] = fileparts(P{n});
    fname         = fullfile(dir,[prefix nam ext]);
    
    Nii = nifti(fname);
    img = squeeze(Nii.dat(:,:,:));
    
    subplot(F1,F2,n);
    imagesc(img); axis image xy off; title(['n=' num2str(n)]); colormap(gray);
end  
drawnow

%% Clean-up
for n=1:N
   delete(P{n});
end