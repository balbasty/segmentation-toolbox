function B = log_spatial_priors(B,wp)
B   = bsxfun(@times,B,wp);
B   = log(bsxfun(@times,B,1./sum(B,2)));
%=======================================================================