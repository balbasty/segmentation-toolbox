function msk = get_msk(f)
msk = isfinite(f) & (f~=0) & (f~=-3024) & (f~=-1500) & (f~=3071) & (f~=max(f(:))) & (f~=min(f(:)));