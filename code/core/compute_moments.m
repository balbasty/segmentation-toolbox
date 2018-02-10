function [mom,ll,mgm] = compute_moments(buf,K,mg,gmm,wp,lkp,K_lab)
    
nz = numel(buf);
N  = numel(buf(1).f);

if nargout==3
    Kb   = numel(wp);
    mgm  = zeros(1,Kb);
end

ll  = 0;
mom = moments_struct(K,N);
for z=1:nz
    if ~buf(z).Nm, continue; end
    
    B = double(buf(z).dat);
    if nargout==3
        s   = 1./(B*wp');
        mgm = mgm + s'*B;
    end
    
    cr                             = NaN(numel(buf(z).msk{1}),N);
    for n=1:N, cr(buf(z).msk{n},n) = double(buf(z).f{n}).*double(buf(z).bf{n}); end    
    
    if nargin>2
        [q,dll] = latent(buf(z).f,buf(z).bf,mg,gmm,B,lkp,wp,buf(z).msk,buf(z).code,K_lab,cr);
        ll      = ll + dll;
    else
        msk1   = buf(z).code>0;
        q      = NaN(numel(buf(z).msk{1}),N);
        for k=1:K
            q(msk1,k) = double(buf(z).dat(:,k));
        end
    end
    
    mom = spm_SuffStats(cr,q,mom,buf(z).code);
end
%=======================================================================