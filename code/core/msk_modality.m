function msk = msk_modality(f,modality,trunc_ct)
if strcmp(modality,'MRI'),    
    msk = isfinite(f) & (f~=0);
elseif strcmp(modality,'CT'), 
    if trunc_ct
        msk = isfinite(f) & (f>0) & (f<=100);  
    else
        msk = isfinite(f) & (f~=min(f(:))) & (f~=0);         
    end    
end
%==========================================================================