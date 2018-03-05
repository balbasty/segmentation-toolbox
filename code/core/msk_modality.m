function msk = msk_modality(f,modality,trunc_ct)
if strcmp(modality,'MRI'),    
    msk = isfinite(f) & (f~=0);
elseif strcmp(modality,'CT'), 
    if ~isempty(trunc_ct)
        msk = isfinite(f) & (f>=trunc_ct(1)) & (f<=trunc_ct(2));  
    else
        msk = isfinite(f) & (f~=min(f(:))) & (f~=0);         
    end    
end
%==========================================================================
