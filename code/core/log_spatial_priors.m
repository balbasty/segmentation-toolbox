function B = log_spatial_priors(B,wp,wp1)
if isempty(wp), wp = ones(1,size(B,2)); end
if nargin<3,    wp1 = 1; end

B = bsxfun(@times,B,wp);
B = log(wp1*bsxfun(@times,B,1./sum(B,2)));
%=======================================================================