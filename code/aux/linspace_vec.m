function y = linspace_vec(x1,x2,n)

if numel(x1)==1 && numel(x2)==1
    y = linspace(x1,x2,n);
    return
end

x1 = squeeze(x1); x2 = squeeze(x2);

if ndims(x1)~= ndims(x2) || any(size(x1)~= size(x2))
    error('d1 and d2 must have the same number of dimension and the same size'),
end

NDim = ndims(x1);
if NDim==2 && any(size(x1)==1)
    NDim = NDim-1;
    if all(size(x1)==1)
        NDim = 0;
    end
end

pp      = (0:n-2)./(floor(n)-1);

Sum1 = kron(x1, ones(1,n-1));
Sum2 = kron((x2-x1), pp);
y    = cat(NDim+1, Sum1  + Sum2, shiftdim(x2, size(x1, 1)==1 ));