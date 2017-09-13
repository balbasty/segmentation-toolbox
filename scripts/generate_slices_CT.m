close all; clc;

addpath('../util/')
addpath('../core/')
addpath('../preproc/')

% Parameters
datadir  = '/home/mbrud/Data/CT-healthy/';
slicedir = '/home/mbrud/Dropbox/PhD/Data/CT-healthy-2D/';

%% Get triplets of same subject images
f = dir(fullfile(datadir,'*.nii'));
N = numel(f);
D = 1;

P = cell(N,1);
for i=1:N
    P{i} = fullfile(datadir,f(i).name);
end

%% Do processing
if (exist(slicedir,'dir') == 0)
    mkdir(slicedir);
end

for i=1:N
    fprintf('%d ',i); 

    f = fullfile(slicedir,['S' num2str(i)]);

    if exist(f,'dir')
        rmdir(f,'s');
    end

    if (exist(f,'dir') == 0)
        mkdir(f);
    end 

    Ptmp = cell(D,1);
    for i2=1:D           
        copyfile(P{i,i2},f);

        [~,nam,ext] = fileparts(P{i,i2});
        Ptmp{i2}    = fullfile(f,[nam ext]);

        nm_reorient(Ptmp{i2},1.5,1);
        reset_origin(Ptmp{i2});
        rigid_align(Ptmp{i2});  
%         atlas_crop(Ptmp{i2}); 
    end  

    % Warp same size
    samp  = 1.5;

    mat  = [];
    dall = [];
    for i2=1:D  
        Nii = nifti(Ptmp{i2});
        mat = cat(3,mat,Nii.mat);
        dm  = size(Nii.dat);
        if numel(dm)==2
           dm(3) = 1;
        end
        dall = cat(1,dall,dm);
    end

    [mat,dall] = compute_avg_mat(mat,dall);

    vx       = sqrt(sum(mat(1:3,1:3).^2));
    st       = samp./vx;
    F        = diag(st);
    F(1:4,4) = 1;
    mat      = mat*F;
    dall      = max(floor(dall./st),1);

    Pw = cell(D,1);
    for i2=1:D  
        Nii  = nifti(Ptmp{i2});
        matf = Nii.mat;        
        img  = Nii.dat(:,:,:);    

        phi = affine_transf(matf\mat,identity(dall));
        img = warp(single(img),single(phi)); 

        [~,nam,ext] = fileparts(Nii.dat.fname);

        Pw{i2} = fullfile(f,['w' nam ext]);    
        create_nii(Pw{i2},img,mat,'float32','wf');
    end

    for i2=1:D
        delete(Ptmp{i2});    
    end

    % Make 2D
    prefix = 'slice';
    for i2=1:D
        V  = spm_vol(Pw{i2});
        dm = V.dim;
        d1 = floor(dm(3)/2) + 10;
        subvol(V,[-inf inf;-inf inf;d1 d1]',prefix);

        [pth,nam,ext] = fileparts(Pw{i2});
        P{i,i2}       = fullfile(pth,[prefix nam ext]);
    end

    for i2=1:D
        delete(Pw{i2});
    end
end
fprintf('\n'); 

%% Visualise
N = min(N,12);

F1  = floor(sqrt(N));
F2  = ceil(N/F1);      
cnt = 1;
for i1=1:N
    for i2=1:D
        Nii = nifti(P{i1,i2});
        img = squeeze(Nii.dat(:,:,:));

        subplot(F1,F2,i1);
        imagesc(img); axis image xy off; title(['n=' num2str(i1)]); colormap(gray); colorbar;
        
        cnt = cnt + 1;
    end
end  
drawnow