function msk = get_msk(f)
msk = isfinite(f(:)) & f(:)~=0 & f(:)~=min(f(:)) & f(:)~=max(f(:));
% msk = ones(size(f(:)),'logical');
%==========================================================================