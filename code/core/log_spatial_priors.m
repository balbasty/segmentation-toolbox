function B = log_spatial_priors(B,wp)
B = bsxfun(@times,B,wp);
B = max(B,eps);
B = log(bsxfun(@times,B,1./sum(B,2)));
%==========================================================================