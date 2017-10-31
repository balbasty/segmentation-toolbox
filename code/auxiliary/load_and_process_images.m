function [V,N,S] = load_and_process_images(imdir,load_from_tmpdir,preproc,runpar,run_on_2D)

if load_from_tmpdir
    % Read already processed image data into a cell array object
    imdir{1} = imdir{5};
    
    if ~run_on_2D
        sanity_check_dir(imdir);        
        [V,N,S] = get_V(imdir);            
    else
        imdir_2D = [imdir{1} '_2D'];                  
        if exist(imdir_2D,'dir')
            imdir{1} = imdir_2D;
            sanity_check_dir(imdir);
            [V,N,S] = get_V(imdir);
        else
            error('exist(imdir_2D,''dir'')')
        end
    end
else
    % Read image data into a cell array object
    sanity_check_dir(imdir);
    
    [V,N,S] = get_V(imdir);   
    
    % Copy input images to temporary directory
    V = copy2tmpdir(V,imdir{5});        
    
    if preproc.create_2D
        % Create a temporary directory to store 2D versions of processed
        % images
        imdir_2D = [imdir{5} '_2D'];
        
        if exist(imdir_2D,'dir')
            rmdir(imdir_2D,'s');
        end
        mkdir(imdir_2D);
    end
       
    % Process the input images in parallel
    parfor (s=1:S,runpar)
        fprintf('s=%d\n',s); 
        
        if N>1
            % If multi-channel data, register and reslice to image with
            % largest volume
            V{s} = reg_and_reslice(V{s});
        end
    
        Affine = [];
        for n=1:N            
            if preproc.realign
                % Reset origin and align to MNI space
                vx = sqrt(sum(V{s}(n).mat(1:3,1:3).^2));
                
                try
                    nm_reorient(V{s}(n).fname,vx,0);
                    
                    [pth,nam,ext] = fileparts(V{s}(n).fname);
                    delete(V{s}(n).fname);

                    nfname  = fullfile(pth,['ro_' nam ext]);
                    V{s}(n) = spm_vol(nfname);
                catch
                    warning('nm_reorient')
                    disp(V{s}(n).fname);
                end                

                try
                    reset_origin(V{s}(n).fname);
                catch
                    warning('reset_origin')
                    disp(V{s}(n).fname);
                end
                
                try
                    rigid_align(V{s}(n).fname);
                catch
                    warning('rigid_align')
                    disp(V{s}(n).fname);
                end
            end
            
            if preproc.crop
                % Remove data outside of the head     
                try
                    if n==1
                        Affine = atlas_crop(V{s}(n).fname,[],'sv_');
                    else
                        atlas_crop(V{s}(n).fname,Affine,'sv_');
                    end

                    [pth,nam,ext] = fileparts(V{s}(n).fname);
                    delete(V{s}(n).fname);

                    nfname  = fullfile(pth,['sv_' nam ext]);
                    V{s}(n) = spm_vol(nfname);
                catch
                    warning('atlas_crop')
                    disp(V{s}(n).fname);
                end 
            end
            
            if preproc.denoise && strcmp(imdir{4},'CT')
                % Denoise using L2-TV (ADMM)
                try
                    denoise_img(V{s}(n).fname);                
                    
                    [pth,nam,ext] = fileparts(V{s}(n).fname);
                    delete(V{s}(n).fname);

                    nfname  = fullfile(pth,['den_' nam ext]);
                    V{s}(n) = spm_vol(nfname);
                catch
                    warning('denoise_img')
                    disp(V{s}(n).fname);
                end                
            end
            
            if preproc.create_2D
                % Save a 2D slice of the processed image
                dm = V{s}(n).dim;
                if dm(3)>1
                    nz = floor(dm(3)/2) + 1; % z-dimension to slice
                    
                    try
                        subvol(V{s}(n),[-inf inf;-inf inf;nz nz]','2D_');
                        
                        [pth,nam,ext] = fileparts(V{s}(n).fname);               
                        nfname        = fullfile(pth,['2D_' nam ext]);
                        
                        reset_origin(nfname);
                        
                        sdir = fullfile(imdir_2D,['S' num2str(s)]);
                        if ~exist(sdir,'dir')
                            mkdir(sdir);
                        end                        

                        movefile(nfname,sdir);
                    catch
                        warning('Create 2D slice')
                        disp(V{s}(n).fname);
                    end  
                end                                
            end
        end
    end
    
    if run_on_2D
        imdir{1} = imdir_2D;
        V        = get_V(imdir);
    end
end

for s=1:S
    for n=1:N
        V{s}(n).brain_is = imdir{2};
        V{s}(n).descrip  = imdir{4};
    end    
end
%==========================================================================

%=======================================================================
function V = copy2tmpdir(V,tpmdir)
S = numel(V);
N = numel(V{1});

if exist(tpmdir,'dir')
    rmdir(tpmdir,'s');
end
mkdir(tpmdir);

nV = cell(1,S);
for s=1:S
    fname = V{s}.fname;
    pth   = fileparts(fname);
    
    splitstr = strsplit(pth,'/');
    
    f = fullfile(tpmdir,splitstr{end});
    if exist(f,'dir')
        rmdir(f,'s');
    end
    mkdir(f);

    copyfile(pth,f);
    
    for n=1:N
        [~,nam,ext] = fileparts(V{s}(n).fname);
        nV{s}(n)     = spm_vol(fullfile(f,[nam ext]));
    end
end
V = nV;
%=======================================================================

%==========================================================================
function [V,N,S] = get_V(imdir)

S_requested = imdir{3};

folder    = dir(imdir{1});  % folder with subfolders containing multi-channel data of subjects
folder    = folder(3:end);
dirflag   = [folder.isdir];

subfolder = folder(dirflag);   % subfolders (S1,S2,...)
S1        = numel(subfolder);
if S_requested>S1
    S = S1;
else
    S = S_requested;
end

files = dir(fullfile(imdir{1},subfolder(1).name,'*.nii'));
N     = numel(files);

V = cell(1,S);
for s=1:S
    folder = fullfile(imdir{1},subfolder(s).name);
    files  = dir(fullfile(folder,'*.nii'));
    for n=1:N
        V{s}(n) = spm_vol(fullfile(folder,files(n).name));
    end        
end 
fprintf('Loading data from %d subjects having %d channels each\n',S,N); 
%==========================================================================

%==========================================================================
function sanity_check_dir(imdir)
if numel(imdir)~=5
    error('numel(imdir)~=5')
end

folder    = dir(imdir{1}); 
folder    = folder(3:end);
dirflag   = [folder.isdir];

subfolder = folder(dirflag);
S         = numel(subfolder);
if S==0
    error('S==0')
end

files = dir(fullfile(imdir{1},subfolder(1).name,'*.nii'));
N0    = numel(files);
for s=2:S
    files = dir(fullfile(imdir{1},subfolder(s).name,'*.nii'));
    N1    = numel(files);
    if N0~=N1
        error('N0~=N1')
    end
end
%==========================================================================