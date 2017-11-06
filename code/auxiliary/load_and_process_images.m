function [V,M,S,N] = load_and_process_images(imobj,preproc,runpar)

M = numel(imobj);
V = cell(1,M);
S = zeros(1,M);
N = zeros(1,M);

if ~preproc.do_preproc
    % Load S images, of N channels, into a cell array object (V)  
    for m=1:M
        [V{m},S(m),N(m)]  = get_V(imobj{m});
    end
else   
    if preproc.is_DICOM
        % imdhir points to a folder structure of DICOMS, which is converted
        % to NIFTIs
%         if N>1
%            error('N>1'); 
%         end
        
        dirNII = './DICOM2NII';
        if exist(dirNII,'dir')
            rmdir(dirNII,'s');
        end
        mkdir(dirNII);
        
        search_and_convert_dcm(imobj{1},dirNII);
        
        niis = dir(fullfile(dirNII,'*.nii'));
        for s=1:numel(niis)
            dirs = fullfile(dirNII,['S' num2str(s)]);            
            mkdir(dirs);
        
            movefile(fullfile(dirNII,niis(s).name),dirs);
        end
        
        imobj{1} = dirNII;
    end
    
    % Read image data into a cell array object        
    [V,S,N] = get_V(imobj);   
    
    % Copy input images to temporary directory
    [V,imdir_2D] = cpy2imdir(V);            
       
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
            
            if preproc.denoise && strcmp(imobj{3},'CT')
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
%==========================================================================

%=======================================================================
function [V,imdir_2D] = cpy2imdir(V)
imdir = './ims';

S = numel(V);
N = numel(V{1});

if exist(imdir,'dir')
    rmdir(imdir,'s');
end
mkdir(imdir);

imdir_2D = './ims_2D';
if exist(imdir_2D,'dir')
    rmdir(imdir_2D,'s');
end
mkdir(imdir_2D);
    
nV = cell(1,S);
for s=1:S
    fname = V{s}.fname;
    odir   = fileparts(fname);
    
    ndir = fullfile(imdir,['S' num2str(s)]);
    if exist(ndir,'dir')
        rmdir(ndir,'s');
    end
    mkdir(ndir);

    copyfile(odir,ndir);
    
    for n=1:N
        [~,nam,ext] = fileparts(V{s}(n).fname);
        nV{s}(n)    = spm_vol(fullfile(ndir,[nam ext]));
    end
end
V = nV;
%=======================================================================