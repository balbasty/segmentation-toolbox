function s = logsumexp(b,dim)
B         = nanmax(b,[],dim);
dims      = ones(1,ndims(b));
dims(dim) = size(b,dim);
b         = b - repmat(B, dims);
s         = B + log(nansum(exp(b),dim));
%=======================================================================