close all; clc;

addpath('../util/')
addpath('../core/')
addpath('../preproc/')

% Parameters
datadir  = '/home/smajjk/Data/IXI-full/';
regdir   = '/home/smajjk/Data/IXI-reg/';
slicedir = '/home/smajjk/Data/IXI-slices/';

reg   = false;
slice = true;

%% Get triplets of same subject images
d = dir(datadir);
d = d(3:end);
D = numel(d);

fnames = cell(D,1);
for i=1:D
    fnames{i} = dir(fullfile(d(i).folder,d(i).name,'*.nii'));
end

P = {};
for i1=1:numel(fnames{1})
    
    tag1 = fnames{1}(i1).name(1:6);

    Ptmp{1} = fullfile(fnames{1}(i1).folder,fnames{1}(i1).name);
    exists = 0;
    for i2=2:D        
        for i3=1:numel(fnames{i2})

            tag2 = fnames{i2}(i3).name(1:6);

            if strcmp(tag1,tag2)
               exists = 1; 
               
               Ptmp{end + 1} = fullfile(fnames{i2}(i3).folder,fnames{i2}(i3).name);
               
               break;
            end
        end
    end
    
    if exists
        P{end + 1} = Ptmp;
    end
    
    Ptmp = {};
end

%% Co-register triplets
if reg
    if (exist(regdir,'dir') == 0)
        mkdir(regdir);
    end

    prefix = 'reg';
    N      = numel(P);
    % Preg   = cell(N,D);
    for i=1:N
        Ptmp = {};

        ref = P{i}(1);

        Ptmp{end + 1} = ref{1};

        for i2=1:numel(P{i}) - 1
            source = P{i}(i2 + 1);

            [pth,nam,ext] = fileparts(source{1});
            Ptmp{end + 1} = fullfile(pth,[prefix nam ext]); 

            matlabbatch = [];
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = ref;
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = source;
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = prefix;

            output_list = spm_jobman('run',matlabbatch);
        end    

        f = fullfile(regdir,['S' num2str(i)]);

        if exist(f,'dir')
            rmdir(f,'s');
        end
        if (exist(f,'dir') == 0)
            mkdir(f);
        end 

        for i2=1:D           
            copyfile(Ptmp{i2},f);

            [~,nam,ext] = fileparts(Ptmp{i2});
            P{i,i2}     = fullfile(f,[nam ext]);
    %         Preg{i,i2}  = fullfile(f,[nam ext]);
        end    

        for i2=1:D-1         
            delete(Ptmp{i2 + 1});
        end 
    end
end

%% Create slices
if slice
    if (exist(slicedir,'dir') == 0)
        mkdir(slicedir);
    end

    NN = N;

    % Pslice = cell(NN,D);
    for i=1:NN
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

            nm_reorient(Ptmp{i2},1,1);
            reset_origin(Ptmp{i2});
            rigid_align(Ptmp{i2});  
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
            img    = Nii.dat(:,:,:);    

            phi = affine_transf(matf\mat,identity(dall));
            img = warp(img,phi); 

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
    %         Pslice{i,i2}  = fullfile(pth,[prefix nam ext]);
        end

        for i2=1:D
            delete(Pw{i2});
        end
    end
    fprintf('\n'); 
end

%% Visualise
NN = min(NN,6);

F1  = floor(sqrt(NN));
F2  = ceil(NN/F1);      
cnt = 1;
for i1=1:NN
    for i2=1:D
        Nii = nifti(Pslice{i1,i2});
        img = squeeze(Nii.dat(:,:,:));

        subplot(NN,D,cnt);
        imagesc(img); axis image xy off; title(['n=' num2str(i1)]); colormap(gray);
        
        cnt = cnt + 1;
    end
end  
drawnow