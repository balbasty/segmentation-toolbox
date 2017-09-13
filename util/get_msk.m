function msk = get_msk(f,ct)
if nargin<2, ct = false; end

if ~ct
    msk = isfinite(f(:)) & f(:)~=0;
else
    msk = isfinite(f(:)) & f(:)~=0 & f(:)>-1100;
end
%==========================================================================