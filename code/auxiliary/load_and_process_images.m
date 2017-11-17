function [V,M,S,N,labels] = load_and_process_images(imobj,preproc,dir_data,num_workers)

M      = numel(imobj);
V      = cell(1,M);
S      = zeros(1,M);
N      = zeros(1,M);
labels = cell(1,M);

if ~preproc.do_preproc
    % Load S images, of N channels, into a cell array object (V)  
    for m=1:M
        [V{m},S(m),N(m),labels{m}]  = get_V(imobj{m});
        
        if preproc.rem_corrupted && strcmp(imobj{m}{3},'CT')
            [V{m},S(m),labels{m}] = rem_corrupted_ims(V{m},labels{m},num_workers);
        end
    end
else   
    for m=1:M
        if preproc.is_DICOM
            % imobj points to a folder structure of DICOMS, which is converted
            % to NIFTIs
            warning('DICOM conversion started, assumes single-channel data..')
            
            dirNII = fullfile(dir_Data,'DICOM2NII');
            if m>1
                dirNII = [dirNII '_m' num2str(m)];
            end
            if exist(dirNII,'dir')
                rmdir(dirNII,'s');
            end
            mkdir(dirNII);

            search_and_convert_dcm(imobj{m}{1},dirNII);

            niis = dir(fullfile(dirNII,'*.nii'));
            for s=1:numel(niis)
                dirs = fullfile(dirNII,['S' num2str(s)]);            
                mkdir(dirs);

                movefile(fullfile(dirNII,niis(s).name),dirs);
            end

            imobj{m}{1} = dirNII;
        end

        % Read image data into a cell array object        
        [V{m},S(m),N(m),labels{m}] = get_V(imobj{m});   

        % Copy input images to temporary directory
        imdir = fullfile(dir_data,'ims');
        if m>1
            imdir = [imdir '_m' num2str(m)];
        end
        [V{m},imdir_2D] = cpy2imdir(V{m},imdir);            
        
        % Process the input images in parallel
        V_m = V{m};
        parfor (s=1:S(m),num_workers)
            fprintf('s=%d\n',s); 

            if N(m)>1
                % If multi-channel data, register and reslice to image with
                % largest volume
                V_m{s} = reg_and_reslice(V_m{s});
            end

            Affine = [];
            for n=1:N(m)            
                if preproc.realign
                    % Reset origin and align to MNI space
                    vx = sqrt(sum(V_m{s}(n).mat(1:3,1:3).^2));

                    try
                        nm_reorient(V_m{s}(n).fname,vx,0);

                        [pth,nam,ext] = fileparts(V_m{s}(n).fname);
                        delete(V_m{s}(n).fname);

                        nfname  = fullfile(pth,['ro_' nam ext]);
                        V_m{s}(n) = spm_vol(nfname);
                    catch
                        warning('nm_reorient')
                        disp(V_m{s}(n).fname);
                    end                

                    try
                        reset_origin(V_m{s}(n).fname);
                    catch
                        warning('reset_origin')
                        disp(V_m{s}(n).fname);
                    end

                    try
                        rigid_align(V_m{s}(n).fname);
                    catch
                        warning('rigid_align')
                        disp(V_m{s}(n).fname);
                    end
                end

                if preproc.crop
                    % Remove data outside of the head     
                    try
                        if n==1
                            Affine = atlas_crop(V_m{s}(n).fname,[],'sv_',preproc.crop_neck);
                        else
                            atlas_crop(V_m{s}(n).fname,Affine,'sv_',preproc.crop_neck);
                        end

                        [pth,nam,ext] = fileparts(V_m{s}(n).fname);
                        delete(V_m{s}(n).fname);

                        nfname  = fullfile(pth,['sv_' nam ext]);
                        V_m{s}(n) = spm_vol(nfname);
                    catch
                        warning('atlas_crop')
                        disp(V_m{s}(n).fname);
                    end 
                end

                if preproc.denoise && strcmp(imobj{m}{3},'CT')
                    % Denoise using L2-TV (ADMM)
                    try
                        denoise_img(V_m{s}(n).fname);                

                        [pth,nam,ext] = fileparts(V_m{s}(n).fname);
                        delete(V_m{s}(n).fname);

                        nfname  = fullfile(pth,['den_' nam ext]);
                        V_m{s}(n) = spm_vol(nfname);
                    catch
                        warning('denoise_img')
                        disp(V_m{s}(n).fname);
                    end                
                end

                % Save a 2D slice of the processed image
                dm = V_m{s}(n).dim;
                if dm(3)>1
                    nz = floor(dm(3)/2) + 1; % z-dimension to slice

                    try
                        subvol(V_m{s}(n),[-inf inf;-inf inf;nz nz]','2D_');

                        [pth,nam,ext] = fileparts(V_m{s}(n).fname);               
                        nfname        = fullfile(pth,['2D_' nam ext]);

                        reset_origin(nfname);

                        sdir = fullfile(imdir_2D,['S' num2str(s)]);
                        if ~exist(sdir,'dir')
                            mkdir(sdir);
                        end                        

                        movefile(nfname,sdir);
                    catch
                        warning('Create 2D slice')
                        disp(V_m{s}(n).fname);
                    end  
                end 
            end
        end
        V{m} = V_m;

        if preproc.rem_corrupted && strcmp(imobj{m}{3},'CT')
            [V{m},S(m),labels{m}] = rem_corrupted_ims(V{m},labels{m},num_workers);
        end
    end
end
%==========================================================================

%=======================================================================
function [V,imdir_2D] = cpy2imdir(V,imdir)
S = numel(V);
N = numel(V{1});

if exist(imdir,'dir')
    rmdir(imdir,'s');
end
mkdir(imdir);

imdir_2D = [imdir '_2D'];
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