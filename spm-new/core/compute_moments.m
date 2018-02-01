function [mom,ll,mgm] = compute_moments(buf,K,mg,gmm,wp,lkp)
    
nz = numel(buf);
N  = numel(buf(1).f);

if nargout==3
    Kb   = numel(wp);
    mgm  = zeros(1,Kb);
end

ll  = 0;
mom = moments_struct(K,N);
for z=1:nz
    if ~buf(z).nm, continue; end
    
    B = double(buf(z).dat);
    if nargout==3
        s   = 1./(B*wp');
        mgm = mgm + s'*B;
    end
    
    cr                  = zeros(numel(buf(z).f{1}),N);
    for n=1:N, cr(:,n)  = double(buf(z).f{n}.*buf(z).bf{n}); end
    
    if nargin>2
        [q,dll] = latent(buf(z).f,buf(z).bf,mg,gmm,B,lkp,wp,buf(z).code,cr);
        ll      = ll + dll;
    else
        q       = double(buf(z).dat);
    end
    
    mom = spm_SuffStats(cr,q,mom,buf(z).code);
end
%=======================================================================