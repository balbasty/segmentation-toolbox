function resp = init_resp(obj,lkp,d)

if ~isfield(obj,'resp')    
    K        = numel(lkp.part);
    obj.resp = struct;
    if ~isfield(obj.resp,'current')
        obj.resp.current = nifti;
        for k=1:K
            fname               = fullfile(obj.dir_seg,['resp-current' num2str(k) '.nii']);    
            obj.resp.current(k) = spm_misc('create_nii',fname,ones(d,'single')/K,eye(4),spm_type('float32'),'resp-current');
        end
    end
    
    if ~isfield(obj.resp,'search')
        obj.resp.search = nifti;
        for k=1:K
            fname               = fullfile(obj.dir_seg,['resp-search' num2str(k) '.nii']);    
            obj.resp.search(k) = spm_misc('create_nii',fname,ones(d,'single')/K,eye(4),spm_type('float32'),'resp-search');
        end
    end        
    
    resp = obj.resp;
else
    resp = obj.resp;
end
%==========================================================================