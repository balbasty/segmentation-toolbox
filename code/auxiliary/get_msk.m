function msk = get_msk(f,descrip)
if strcmp(descrip,'CT')
    msk = isfinite(f) & (f~=0) & (f~=-3024) & (f~=-1500) & (f~=3071) & (f~=max(f(:))) & (f~=min(f(:)));
elseif strcmp(descrip,'MRI')
    msk = isfinite(f) & (f~=0); 
end