function msk = msk_modality(f,modality)
if strcmp(modality,'MRI'),    
    msk = isfinite(f) & (f~=0);
elseif strcmp(modality,'CT'), 
    msk = isfinite(f) & (f~=min(f(:))) & (f~=0);            
end
%==========================================================================