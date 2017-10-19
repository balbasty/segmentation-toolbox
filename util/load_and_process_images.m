function [V,N,S] = load_and_process_images(imdir,tmpdir,load_from_tmpdir,preproc)

if load_from_tmpdir
    % Read already processed image data into a cell array object
    [V,N,S] = get_V(tmpdir);
else
    % Read image data into a cell array object
    [V,N,S] = get_V(imdir);

    % Copy input images to temporary directory
    V = copy2tmpdir(V,tmpdir{1}{1});        
    
    if preproc.create_2D
        % Create a temporary directory to store 2D versions of processed
        % images
        dir       = fileparts(tmpdir{1}{1});        
        dirTPM_2D = [dir '_2D'];
        
        if exist(dirTPM_2D,'dir')
            rmdir(dirTPM_2D,'s');
        end
        mkdir(dirTPM_2D);
    end
    
    if N>1
        % If multi-channel data, register and reslice to same size
        V = reg_and_reslice(V);
    end
    
    % Process the input images in parallel
    parfor s=1:S   
        for n=1:N            
            if preproc.realign
                % Reset origin and align to MNI space
                vx = sqrt(sum(V{s}(n).mat(1:3,1:3).^2));
                nm_reorient(V{s}(n).fname,vx,0);

                [pth,nam,ext] = fileparts(V{s}(n).fname);
                delete(V{s}(n).fname);

                nfname  = fullfile(pth,['ro_' nam ext]);
                V{s}(n) = spm_vol(nfname);

                reset_origin(V{s}(n).fname);
                
                rigid_align(V{s}(n).fname);
            end
            
            if preproc.crop
                % Remove data outside of the head
                if N>1
                    warning('preproc.crop need to be fixed for N>1')                   
                else                
                    atlas_crop(V{s}(n).fname,[],'sv_');

                    [pth,nam,ext] = fileparts(V{s}(n).fname);
                    delete(V{s}(n).fname);

                    nfname  = fullfile(pth,['sv_' nam ext]);
                    V{s}(n) = spm_vol(nfname);
                end
            end
            
            if preproc.denoise
                % Denoise using L2-TV (ADMM)
                denoise_img(V{s}(n).fname);                

                [pth,nam,ext] = fileparts(V{s}(n).fname);
                delete(V{s}(n).fname);

                nfname  = fullfile(pth,['den_' nam ext]);
                V{s}(n) = spm_vol(nfname);
            end
            
            if preproc.create_2D
                % Save a 2D slice of the processed image
                dm = V{s}(n).dim;
                if dm(3)>1
                    d1 = floor(dm(3)/2) + 1;
                    subvol(V{s}(n),[-inf inf;-inf inf;d1 d1]','2D_');
                    
                    [pth,nam,ext] = fileparts(V{s}(n).fname);               
                    nfname        = fullfile(pth,['2D_' nam ext]);

                    reset_origin(nfname);
                    
                    sdir = fullfile(dirTPM_2D,['S' num2str(s)]);
                    if exist(sdir,'dir')
                        rmdir(sdir,'s');
                    end
                    mkdir(sdir);
        
                    movefile(nfname,sdir);
                end                                
            end
        end
    end
end

