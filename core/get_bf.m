function bf = get_bf(chan,C,d)
bf  = zeros(d,'single');
for c=1:C
    for z=1:d(3)
        bfz       = transf(chan(c).B1,chan(c).B2,chan(c).B3(z,:),chan(c).T);        
        bf(:,:,z) = single(exp(bfz));
    end
end
bf                        = bf(:);
bf(~isfinite(bf) | bf==0) = 1;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
%=======================================================================